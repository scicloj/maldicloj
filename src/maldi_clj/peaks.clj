(ns maldi-clj.peaks
  "High-performance peak detection algorithms for mass spectrometry data"
  (:require [maldi-clj.spectrum :as spectrum]
            [maldi-clj.baseline :as baseline]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [malli.core :as m]))

;; ==============================================================================
;; PEAK DATA STRUCTURES
;; ==============================================================================

;; Schema for individual peak
(def peak-schema
  [:map
   [:mz number?]
   [:intensity number?]
   [:snr {:optional true} number?]
   [:noise {:optional true} number?]
   [:index {:optional true} int?]])

;; Schema for peak collection (MassPeaks equivalent)
(def mass-peaks-schema
  [:map
   [:id string?]
   [:peaks [:vector peak-schema]]
   [:metadata map?]])

;; ==============================================================================
;; NOISE ESTIMATION ALGORITHMS
;; ==============================================================================

(defn estimate-noise-mad
  "Estimate noise using Median Absolute Deviation (MAD).
   
   MAD is a robust noise estimator:
   MAD = median(|x_i - median(x)|)
   noise = 1.4826 * MAD  (for normal distribution assumption)
   
   Parameters:
   - intensities: intensity values (dtype-next container or vector)
   
   Returns estimated noise level"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Efficient dtype-next implementation
    (let [median-val (dfn/median intensities)
          deviations (dfn/abs (dfn/- intensities median-val))
          mad (dfn/median deviations)]
      (* 1.4826 mad))
    ;; Fallback for legacy vectors
    (let [sorted-intensities (sort intensities)
          n (count sorted-intensities)
          median-val (if (even? n)
                       (/ (+ (nth sorted-intensities (/ n 2))
                             (nth sorted-intensities (dec (/ n 2)))) 2.0)
                       (nth sorted-intensities (/ n 2)))
          deviations (map #(Math/abs (- % median-val)) intensities)
          sorted-deviations (sort deviations)
          mad (if (even? n)
                (/ (+ (nth sorted-deviations (/ n 2))
                      (nth sorted-deviations (dec (/ n 2)))) 2.0)
                (nth sorted-deviations (/ n 2)))]
      (* 1.4826 mad))))

(defn estimate-noise-local-mad
  "Estimate local noise using sliding window MAD.
   
   For each point, estimate noise in a local window around it.
   This handles varying noise levels across the spectrum.
   
   Parameters:
   - intensities: intensity values
   - window-half-size: half-size of the sliding window (default 25)
   
   Returns vector of local noise estimates"
  [intensities & {:keys [window-half-size] :or {window-half-size 25}}]
  (if (dtype/reader? intensities)
    ;; Efficient dtype-next implementation
    (let [n (dtype/ecount intensities)
          result (dtype/make-container :jvm-heap :float32 n)]
      (dotimes [i n]
        (let [start (max 0 (- i window-half-size))
              end (min n (+ i window-half-size 1))
              window (dtype/sub-buffer intensities start (- end start))
              noise (estimate-noise-mad window)]
          (dtype/set-value! result i noise)))
      result)
    ;; Fallback for legacy vectors
    (let [n (count intensities)]
      (mapv (fn [i]
              (let [start (max 0 (- i window-half-size))
                    end (min n (+ i window-half-size 1))
                    window (subvec (vec intensities) start end)]
                (estimate-noise-mad window)))
            (range n)))))

;; ==============================================================================
;; LOCAL MAXIMA DETECTION
;; ==============================================================================

(defn find-local-maxima
  "Find local maxima in intensity data.
   
   A point is a local maximum if it's greater than or equal to all points
   in its neighborhood window.
   
   Parameters:
   - intensities: intensity values
   - window-half-size: half-size of the window for local comparison (default 1)
   - min-intensity: minimum intensity threshold (default 0)
   
   Returns vector of indices where local maxima occur"
  [intensities & {:keys [window-half-size min-intensity]
                  :or {window-half-size 1 min-intensity 0}}]
  (if (dtype/reader? intensities)
    ;; Efficient dtype-next implementation
    (let [n (dtype/ecount intensities)
          maxima-indices (transient [])]
      (dotimes [i n]
        (let [current-intensity (dtype/get-value intensities i)]
          (when (>= current-intensity min-intensity)
            (let [start (max 0 (- i window-half-size))
                  end (min n (+ i window-half-size 1))
                  is-maximum? (loop [j start
                                     is-max? true]
                                (cond
                                  (>= j end) is-max?
                                  (= j i) (recur (inc j) is-max?)
                                  (> (dtype/get-value intensities j) current-intensity) false
                                  :else (recur (inc j) is-max?)))]
              (when is-maximum?
                (conj! maxima-indices i))))))
      (persistent! maxima-indices))
    ;; Fallback for legacy vectors
    (let [n (count intensities)]
      (filterv (fn [i]
                 (let [current-intensity (nth intensities i)]
                   (when (>= current-intensity min-intensity)
                     (let [start (max 0 (- i window-half-size))
                           end (min n (+ i window-half-size 1))
                           window (subvec (vec intensities) start end)
                           max-in-window (apply max window)]
                       (= current-intensity max-in-window)))))
               (range n)))))

;; ==============================================================================
;; SNR CALCULATION
;; ==============================================================================

(defn calculate-snr
  "Calculate Signal-to-Noise Ratio.
   
   SNR = (signal - noise) / noise
   
   Parameters:
   - intensity: signal intensity
   - noise: estimated noise level
   
   Returns SNR value as double"
  [intensity noise]
  (if (> noise 0)
    (double (/ (- intensity noise) noise))
    Double/POSITIVE_INFINITY))

(defn calculate-snr-vector
  "Calculate SNR for a vector of intensities given noise estimates.
   
   Parameters:
   - intensities: intensity values
   - noise-estimates: noise estimates (single value or vector)
   
   Returns vector of SNR values"
  [intensities noise-estimates]
  (if (dtype/reader? intensities)
    ;; Efficient dtype-next implementation
    (if (number? noise-estimates)
      ;; Single noise value for all points
      (dfn/max (dfn// (dfn/- intensities noise-estimates) noise-estimates) 0.0)
      ;; Vector of noise estimates
      (dfn/max (dfn// (dfn/- intensities noise-estimates) noise-estimates) 0.0))
    ;; Fallback for legacy vectors
    (if (number? noise-estimates)
      (mapv #(max (calculate-snr % noise-estimates) 0.0) intensities)
      (mapv #(max (calculate-snr %1 %2) 0.0) intensities noise-estimates))))

;; ==============================================================================
;; MAIN PEAK DETECTION FUNCTION
;; ==============================================================================

(defn detect-peaks
  "Detect peaks in a mass spectrum using SNR-based filtering.
   
   This is the main peak detection function equivalent to MALDIquant's detectPeaks.
   
   Algorithm:
   1. Apply baseline correction (if not already done)
   2. Estimate noise levels
   3. Find local maxima
   4. Calculate SNR for each maximum
   5. Filter peaks by SNR threshold
   6. Return peak collection
   
   Parameters:
   - spectrum: mass spectrum data structure
   - Options:
     :snr-threshold - minimum SNR for peak detection (default 3.0)
     :window-half-size - half-size for local maxima detection (default 10)
     :noise-method - :mad or :local-mad (default :mad)
     :noise-window - window size for local noise estimation (default 50)
     :min-intensity - minimum intensity threshold (default 0)
     :baseline-corrected? - skip baseline correction if true (default false)
   
   Returns MassPeaks data structure"
  [spectrum & {:keys [snr-threshold window-half-size noise-method noise-window
                      min-intensity baseline-corrected?]
               :or {snr-threshold 3.0
                    window-half-size 10
                    noise-method :mad
                    noise-window 50
                    min-intensity 0
                    baseline-corrected? false}}]

  ;; Input validation
  (when-not (spectrum/valid-spectrum? spectrum)
    (throw (ex-info "Invalid spectrum" {:spectrum spectrum})))

  (when (= 0 (spectrum/spectrum-length spectrum))
    (throw (ex-info "Cannot detect peaks in empty spectrum" {:spectrum spectrum})))

  (let [;; Step 1: Baseline correction (if needed)
        corrected-spectrum (if baseline-corrected?
                             spectrum
                             (baseline/remove-baseline spectrum :snip))

        mz-values (:mz-values corrected-spectrum)
        intensities (:intensities corrected-spectrum)

        ;; Step 2: Noise estimation
        noise-estimates (case noise-method
                          :mad (estimate-noise-mad intensities)
                          :local-mad (estimate-noise-local-mad intensities
                                                               :window-half-size (/ noise-window 2))
                          (throw (ex-info "Unknown noise method" {:method noise-method})))

        ;; Step 3: Find local maxima
        maxima-indices (find-local-maxima intensities
                                          :window-half-size window-half-size
                                          :min-intensity min-intensity)

        ;; Step 4: Calculate SNR and filter
        peaks (reduce
               (fn [acc idx]
                 (let [mz (if (dtype/reader? mz-values)
                            (dtype/get-value mz-values idx)
                            (nth mz-values idx))
                       intensity (if (dtype/reader? intensities)
                                   (dtype/get-value intensities idx)
                                   (nth intensities idx))
                       noise (if (number? noise-estimates)
                               noise-estimates
                               (if (dtype/reader? noise-estimates)
                                 (dtype/get-value noise-estimates idx)
                                 (nth noise-estimates idx)))
                       snr (calculate-snr intensity noise)]

                   ;; Filter by SNR threshold
                   (if (>= snr snr-threshold)
                     (conj acc {:mz mz
                                :intensity intensity
                                :snr snr
                                :noise noise
                                :index idx})
                     acc)))
               []
               maxima-indices)]

    ;; Return MassPeaks structure
    {:id (:id spectrum)
     :peaks (vec peaks)
     :metadata (merge (:metadata spectrum)
                      {:peak-detection {:snr-threshold snr-threshold
                                        :window-half-size window-half-size
                                        :noise-method noise-method
                                        :baseline-corrected? baseline-corrected?
                                        :total-peaks (count peaks)
                                        :total-maxima (count maxima-indices)}})}))

;; ==============================================================================
;; PEAK UTILITIES
;; ==============================================================================

(defn peak-count
  "Get number of peaks in a MassPeaks collection"
  [mass-peaks]
  (count (:peaks mass-peaks)))

(defn get-peak-intensities
  "Extract peak intensities as a vector"
  [mass-peaks]
  (mapv :intensity (:peaks mass-peaks)))

(defn get-peak-mz-values
  "Extract peak m/z values as a vector"
  [mass-peaks]
  (mapv :mz (:peaks mass-peaks)))

(defn get-peak-snr-values
  "Extract peak SNR values as a vector"
  [mass-peaks]
  (mapv :snr (:peaks mass-peaks)))

(defn filter-peaks-by-snr
  "Filter peaks by SNR threshold"
  [mass-peaks min-snr]
  (update mass-peaks :peaks
          (fn [peaks]
            (filterv #(>= (:snr %) min-snr) peaks))))

(defn filter-peaks-by-intensity
  "Filter peaks by intensity threshold"
  [mass-peaks min-intensity]
  (update mass-peaks :peaks
          (fn [peaks]
            (filterv #(>= (:intensity %) min-intensity) peaks))))

(defn peak-statistics
  "Calculate statistics for a peak collection"
  [mass-peaks]
  (let [peaks (:peaks mass-peaks)
        intensities (get-peak-intensities mass-peaks)
        snr-values (get-peak-snr-values mass-peaks)]
    (when (seq peaks)
      {:count (count peaks)
       :intensity-stats {:min (apply min intensities)
                         :max (apply max intensities)
                         :mean (/ (reduce + intensities) (count intensities))}
       :snr-stats {:min (apply min snr-values)
                   :max (apply max snr-values)
                   :mean (/ (reduce + snr-values) (count snr-values))}})))

;; ==============================================================================
;; CONVENIENCE FUNCTIONS
;; ==============================================================================

(defn estimate-noise
  "Convenience function for noise estimation with multiple methods.
   
   Parameters:
   - spectrum-or-intensities: full spectrum or just intensity values
   - method: :mad or :local-mad
   - options: method-specific options
   
   Returns noise estimate(s)"
  [spectrum-or-intensities method & options]
  (let [intensities (if (map? spectrum-or-intensities)
                      (:intensities spectrum-or-intensities)
                      spectrum-or-intensities)]
    (case method
      :mad (estimate-noise-mad intensities)
      :local-mad (apply estimate-noise-local-mad intensities options)
      (throw (ex-info "Unknown noise estimation method" {:method method})))))
