(ns maldi-clj.alignment
  "Peak alignment and binning algorithms for multi-spectrum analysis"
  (:require [maldi-clj.spectrum :as spectrum]
            [maldi-clj.peaks :as peaks]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]

            [malli.core :as m]))

;; ==============================================================================
;; DATA STRUCTURES FOR ALIGNED PEAKS
;; ==============================================================================

;; Schema for binned peaks across multiple spectra
(def binned-peaks-schema
  [:map
   [:bin-id string?]
   [:mz-center number?] ; Consensus m/z for this bin
   [:mz-range [:tuple number? number?]] ; [min-mz max-mz] for this bin
   [:peaks [:vector [:map ; Peak list with spectrum info
                     [:spectrum-id string?]
                     [:mz number?]
                     [:intensity number?]
                     [:snr {:optional true} number?]
                     [:noise {:optional true} number?]]]]
   [:frequency int?] ; Number of spectra containing this peak
   [:metadata map?]])

;; Schema for collection of binned peaks
(def peak-bins-schema
  [:map
   [:bins [:vector binned-peaks-schema]]
   [:total-spectra int?]
   [:binning-parameters map?]
   [:metadata map?]])

;; ==============================================================================
;; M/Z TOLERANCE AND MATCHING UTILITIES
;; ==============================================================================

(defn calculate-mz-tolerance
  "Calculate m/z tolerance based on mode and parameters.
   
   Parameters:
   - mz: m/z value
   - tolerance: tolerance value
   - tolerance-mode: :ppm or :absolute
   
   Returns [min-mz max-mz] tolerance range"
  [mz tolerance tolerance-mode]
  (case tolerance-mode
    :ppm (let [delta (* mz tolerance 1e-6)]
           [(- mz delta) (+ mz delta)])
    :absolute [(- mz tolerance) (+ mz tolerance)]
    (throw (ex-info "Invalid tolerance mode" {:mode tolerance-mode}))))

(defn mz-within-tolerance?
  "Check if two m/z values are within tolerance"
  [mz1 mz2 tolerance tolerance-mode]
  (let [[min-mz max-mz] (calculate-mz-tolerance mz1 tolerance tolerance-mode)]
    (<= min-mz mz2 max-mz)))

(defn find-matching-peaks
  "Find peaks from peak-list that match target-mz within tolerance.
   
   Returns vector of matching peak indices"
  [peak-list target-mz tolerance tolerance-mode]
  (let [mz-values (if (map? (first peak-list))
                    (mapv :mz peak-list)
                    peak-list)]
    (keep-indexed
     (fn [idx mz]
       (when (mz-within-tolerance? target-mz mz tolerance tolerance-mode)
         idx))
     mz-values)))

;; ==============================================================================
;; PEAK BINNING ALGORITHMS
;; ==============================================================================

(defn bin-peaks-strict
  "Bin peaks using strict mode - peaks must match exactly across spectra.
   
   This is the most conservative approach, only creating bins for peaks
   that appear at very similar m/z values across multiple spectra."
  [mass-peaks-collection tolerance tolerance-mode min-frequency]
  (let [;; Collect all unique m/z values across all spectra
        all-peaks (mapcat (fn [mass-peaks]
                            (map #(assoc % :spectrum-id (:id mass-peaks))
                                 (:peaks mass-peaks)))
                          mass-peaks-collection)

        ;; Group peaks by m/z tolerance
        peak-groups (reduce
                     (fn [groups peak]
                       (let [matching-group-idx (first
                                                 (keep-indexed
                                                  (fn [idx group]
                                                    (when (some #(mz-within-tolerance?
                                                                  (:mz peak) (:mz %)
                                                                  tolerance tolerance-mode)
                                                                group)
                                                      idx))
                                                  groups))]
                         (if matching-group-idx
                           (update groups matching-group-idx conj peak)
                           (conj groups [peak]))))
                     []
                     all-peaks)

        ;; Filter groups by minimum frequency and create bins
        valid-bins (keep-indexed
                    (fn [idx group]
                      (let [unique-spectra (set (map :spectrum-id group))
                            frequency (count unique-spectra)]
                        (when (>= frequency min-frequency)
                          (let [mz-values (map :mz group)
                                mz-center (/ (reduce + mz-values) (count mz-values))
                                mz-min (apply min mz-values)
                                mz-max (apply max mz-values)]
                            {:bin-id (str "bin-" idx)
                             :mz-center mz-center
                             :mz-range [mz-min mz-max]
                             :peaks (vec group)
                             :frequency frequency
                             :metadata {:binning-mode :strict}}))))
                    peak-groups)]

    {:bins (vec valid-bins)
     :total-spectra (count mass-peaks-collection)
     :binning-parameters {:tolerance tolerance
                          :tolerance-mode tolerance-mode
                          :min-frequency min-frequency
                          :mode :strict}
     :metadata {:created-at (java.time.Instant/now)}}))

(defn bin-peaks-relaxed
  "Bin peaks using relaxed mode - more permissive m/z matching.
   
   This approach creates bins more liberally, allowing for slight m/z
   variations that might occur due to calibration differences."
  [mass-peaks-collection tolerance tolerance-mode min-frequency]
  (let [;; Create a larger tolerance for relaxed mode
        relaxed-tolerance (* tolerance 2.0)

        ;; Get all peaks with spectrum info
        all-peaks (mapcat (fn [mass-peaks]
                            (map #(assoc % :spectrum-id (:id mass-peaks))
                                 (:peaks mass-peaks)))
                          mass-peaks-collection)

        ;; Sort peaks by m/z for efficient processing
        sorted-peaks (sort-by :mz all-peaks)

        ;; Create bins using sliding window approach
        bins (loop [peaks sorted-peaks
                    current-bin []
                    all-bins []
                    bin-counter 0]
               (if (empty? peaks)
                 ;; Process final bin
                 (if (>= (count (set (map :spectrum-id current-bin))) min-frequency)
                   (conj all-bins current-bin)
                   all-bins)

                 (let [current-peak (first peaks)
                       remaining-peaks (rest peaks)]

                   (if (empty? current-bin)
                     ;; Start new bin
                     (recur remaining-peaks [current-peak] all-bins bin-counter)

                     ;; Check if current peak fits in current bin
                     (let [bin-center-mz (/ (reduce + (map :mz current-bin)) (count current-bin))]
                       (if (mz-within-tolerance? bin-center-mz (:mz current-peak)
                                                 relaxed-tolerance tolerance-mode)
                         ;; Add to current bin
                         (recur remaining-peaks (conj current-bin current-peak) all-bins bin-counter)

                         ;; Start new bin
                         (let [unique-spectra (set (map :spectrum-id current-bin))
                               new-all-bins (if (>= (count unique-spectra) min-frequency)
                                              (conj all-bins current-bin)
                                              all-bins)]
                           (recur remaining-peaks [current-peak] new-all-bins (inc bin-counter)))))))))

        ;; Convert bins to proper format
        formatted-bins (map-indexed
                        (fn [idx bin]
                          (let [mz-values (map :mz bin)
                                mz-center (/ (reduce + mz-values) (count mz-values))
                                mz-min (apply min mz-values)
                                mz-max (apply max mz-values)
                                frequency (count (set (map :spectrum-id bin)))]
                            {:bin-id (str "relaxed-bin-" idx)
                             :mz-center mz-center
                             :mz-range [mz-min mz-max]
                             :peaks (vec bin)
                             :frequency frequency
                             :metadata {:binning-mode :relaxed}}))
                        bins)]

    {:bins (vec formatted-bins)
     :total-spectra (count mass-peaks-collection)
     :binning-parameters {:tolerance tolerance
                          :tolerance-mode tolerance-mode
                          :min-frequency min-frequency
                          :mode :relaxed}
     :metadata {:created-at (java.time.Instant/now)}}))

(defn bin-peaks-reference
  "Bin peaks using reference mode - align all peaks to a reference spectrum.
   
   This approach uses one spectrum as a reference and aligns all other
   peaks to the reference peaks within tolerance."
  [mass-peaks-collection reference-id tolerance tolerance-mode min-frequency]
  (let [;; Find reference spectrum
        reference-spectrum (first (filter #(= (:id %) reference-id) mass-peaks-collection))
        other-spectra (filter #(not= (:id %) reference-id) mass-peaks-collection)]

    (when-not reference-spectrum
      (throw (ex-info "Reference spectrum not found" {:reference-id reference-id})))

    (let [reference-peaks (:peaks reference-spectrum)

          ;; For each reference peak, find matching peaks in other spectra
          bins (keep-indexed
                (fn [idx ref-peak]
                  (let [matching-peaks (reduce
                                        (fn [acc mass-peaks]
                                          (let [matches (filter
                                                         #(mz-within-tolerance?
                                                           (:mz ref-peak) (:mz %)
                                                           tolerance tolerance-mode)
                                                         (:peaks mass-peaks))]
                                            (concat acc
                                                    (map #(assoc % :spectrum-id (:id mass-peaks))
                                                         matches))))
                                        [(assoc ref-peak :spectrum-id reference-id)]
                                        other-spectra)

                        unique-spectra (set (map :spectrum-id matching-peaks))
                        frequency (count unique-spectra)]

                    (when (>= frequency min-frequency)
                      (let [mz-values (map :mz matching-peaks)
                            mz-min (apply min mz-values)
                            mz-max (apply max mz-values)]
                        {:bin-id (str "ref-bin-" idx)
                         :mz-center (:mz ref-peak) ; Use reference m/z as center
                         :mz-range [mz-min mz-max]
                         :peaks (vec matching-peaks)
                         :frequency frequency
                         :metadata {:binning-mode :reference
                                    :reference-spectrum reference-id}}))))
                reference-peaks)]

      {:bins (vec bins)
       :total-spectra (count mass-peaks-collection)
       :binning-parameters {:tolerance tolerance
                            :tolerance-mode tolerance-mode
                            :min-frequency min-frequency
                            :mode :reference
                            :reference-spectrum reference-id}
       :metadata {:created-at (java.time.Instant/now)}})))

;; ==============================================================================
;; MAIN BINNING FUNCTION
;; ==============================================================================

(defn bin-peaks
  "Bin peaks across multiple spectra using specified mode.
   
   This is the main peak binning function equivalent to MALDIquant's binPeaks.
   
   Parameters:
   - mass-peaks-collection: vector of MassPeaks from multiple spectra
   - Options:
     :tolerance - m/z tolerance (default 0.002)
     :tolerance-mode - :ppm or :absolute (default :absolute)
     :method - :strict, :relaxed, or :reference (default :strict)
     :reference-spectrum - spectrum ID for reference mode
     :min-frequency - minimum number of spectra for a valid bin (default 0.5 * total-spectra)
   
   Returns PeakBins structure with aligned peaks"
  [mass-peaks-collection & {:keys [tolerance tolerance-mode method reference-spectrum min-frequency]
                            :or {tolerance 0.002
                                 tolerance-mode :absolute
                                 method :strict
                                 min-frequency nil}}]

  ;; Input validation
  (when (empty? mass-peaks-collection)
    (throw (ex-info "Empty mass peaks collection" {})))

  (let [total-spectra (count mass-peaks-collection)
        actual-min-frequency (or min-frequency (max 1 (int (* 0.5 total-spectra))))]

    (case method
      :strict (bin-peaks-strict mass-peaks-collection tolerance tolerance-mode actual-min-frequency)
      :relaxed (bin-peaks-relaxed mass-peaks-collection tolerance tolerance-mode actual-min-frequency)
      :reference (bin-peaks-reference mass-peaks-collection reference-spectrum
                                      tolerance tolerance-mode actual-min-frequency)
      (throw (ex-info "Invalid binning method" {:method method})))))

;; ==============================================================================
;; PEAK FILTERING FUNCTIONS
;; ==============================================================================

(defn filter-peaks-by-frequency
  "Filter peaks based on frequency of occurrence across spectra.
   
   Parameters:
   - peak-bins: result from bin-peaks
   - min-frequency: minimum frequency threshold
   - frequency-mode: :absolute (count) or :relative (proportion)
   
   Returns filtered peak bins"
  [peak-bins min-frequency & {:keys [frequency-mode] :or {frequency-mode :absolute}}]
  (let [total-spectra (:total-spectra peak-bins)
        threshold (case frequency-mode
                    :absolute min-frequency
                    :relative (int (* min-frequency total-spectra))
                    min-frequency)

        filtered-bins (filter #(>= (:frequency %) threshold) (:bins peak-bins))]

    (assoc peak-bins
           :bins (vec filtered-bins)
           :metadata (assoc (:metadata peak-bins)
                            :filtered-by-frequency true
                            :frequency-threshold threshold
                            :frequency-mode frequency-mode))))

(defn filter-peaks-by-intensity
  "Filter peak bins based on intensity criteria.
   
   Parameters:
   - peak-bins: result from bin-peaks
   - min-intensity: minimum intensity threshold
   - intensity-mode: :any (at least one peak), :all (all peaks), :mean (mean intensity)
   
   Returns filtered peak bins"
  [peak-bins min-intensity & {:keys [intensity-mode] :or {intensity-mode :any}}]
  (let [filtered-bins (filter
                       (fn [bin]
                         (let [intensities (map :intensity (:peaks bin))]
                           (case intensity-mode
                             :any (some #(>= % min-intensity) intensities)
                             :all (every? #(>= % min-intensity) intensities)
                             :mean (>= (/ (reduce + intensities) (count intensities)) min-intensity))))
                       (:bins peak-bins))]

    (assoc peak-bins
           :bins (vec filtered-bins)
           :metadata (assoc (:metadata peak-bins)
                            :filtered-by-intensity true
                            :intensity-threshold min-intensity
                            :intensity-mode intensity-mode))))

;; ==============================================================================
;; PEAK MERGING AND COMBINATION
;; ==============================================================================

(defn merge-mass-peaks
  "Merge multiple MassPeaks collections into a single collection.
   
   Parameters:
   - mass-peaks-collections: vector of MassPeaks
   - merge-strategy: :union (combine all), :intersection (common only)
   
   Returns merged MassPeaks collection"
  [mass-peaks-collections & {:keys [merge-strategy] :or {merge-strategy :union}}]
  (case merge-strategy
    :union (let [all-peaks (mapcat :peaks mass-peaks-collections)
                 all-ids (map :id mass-peaks-collections)
                 merged-metadata (reduce merge (map :metadata mass-peaks-collections))]
             {:id (str "merged-" (clojure.string/join "-" all-ids))
              :peaks (vec all-peaks)
              :metadata (assoc merged-metadata
                               :merge-strategy :union
                               :source-spectra all-ids)})

    :intersection (throw (ex-info "Intersection merge not yet implemented" {}))))

;; ==============================================================================
;; INTENSITY MATRIX GENERATION
;; ==============================================================================

(defn create-intensity-matrix
  "Create intensity matrix from binned peaks.
   
   Rows represent spectra, columns represent m/z bins.
   Missing values are handled according to missing-value-strategy.
   
   Parameters:
   - peak-bins: result from bin-peaks
   - missing-value-strategy: :zero, :na, or :interpolate
   
   Returns map with matrix data structure"
  [peak-bins & {:keys [missing-value-strategy] :or {missing-value-strategy :zero}}]
  (let [bins (:bins peak-bins)
        spectrum-ids (distinct (mapcat #(map :spectrum-id (:peaks %)) bins))

        ;; Create column names from bin centers
        column-names (mapv #(format "mz_%.4f" (:mz-center %)) bins)

        ;; Create matrix data
        matrix-rows (mapv
                     (fn [spectrum-id]
                       (let [row-data (mapv
                                       (fn [bin]
                                         (let [spectrum-peaks (filter #(= (:spectrum-id %) spectrum-id) (:peaks bin))]
                                           (if (empty? spectrum-peaks)
                                             (case missing-value-strategy
                                               :zero 0.0
                                               :na ##NaN
                                               :interpolate 0.0) ; TODO: implement interpolation
                                             ;; Use max intensity if multiple peaks in bin
                                             (apply max (map :intensity spectrum-peaks)))))
                                       bins)]
                         {:spectrum-id spectrum-id
                          :intensities row-data}))
                     spectrum-ids)]

    {:spectrum-ids (vec spectrum-ids)
     :column-names column-names
     :matrix-data matrix-rows
     :dimensions [(count spectrum-ids) (count bins)]
     :metadata {:missing-value-strategy missing-value-strategy
                :created-from-bins (count bins)
                :total-spectra (count spectrum-ids)}}))

;; ==============================================================================
;; UTILITY FUNCTIONS
;; ==============================================================================

(defn bin-statistics
  "Calculate statistics for peak bins"
  [peak-bins]
  (let [bins (:bins peak-bins)
        frequencies (map :frequency bins)
        mz-centers (map :mz-center bins)
        total-peaks (reduce + (map #(count (:peaks %)) bins))]

    {:total-bins (count bins)
     :total-peaks total-peaks
     :frequency-stats {:min (if (seq frequencies) (apply min frequencies) 0)
                       :max (if (seq frequencies) (apply max frequencies) 0)
                       :mean (if (seq frequencies) (/ (reduce + frequencies) (count frequencies)) 0)}
     :mz-range {:min (if (seq mz-centers) (apply min mz-centers) 0)
                :max (if (seq mz-centers) (apply max mz-centers) 0)}
     :coverage (if (pos? (:total-spectra peak-bins))
                 (/ (count bins) (:total-spectra peak-bins))
                 0)}))

(defn get-bin-by-mz
  "Find bin containing the specified m/z value"
  [peak-bins target-mz tolerance tolerance-mode]
  (first (filter #(mz-within-tolerance? (:mz-center %) target-mz tolerance tolerance-mode)
                 (:bins peak-bins))))

(defn export-bins-summary
  "Export summary of binned peaks for analysis"
  [peak-bins]
  (let [bins (:bins peak-bins)]
    (map (fn [bin]
           {:bin-id (:bin-id bin)
            :mz-center (:mz-center bin)
            :mz-range (:mz-range bin)
            :frequency (:frequency bin)
            :total-intensity (reduce + (map :intensity (:peaks bin)))
            :mean-intensity (/ (reduce + (map :intensity (:peaks bin))) (count (:peaks bin)))
            :spectra-present (distinct (map :spectrum-id (:peaks bin)))})
         bins)))
