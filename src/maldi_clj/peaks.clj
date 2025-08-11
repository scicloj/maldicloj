(ns maldi-clj.peaks
  "Peak detection algorithms for MALDI mass spectrometry data.
   Includes optimized implementations for performance-critical operations."
  (:require [clojure.math :as math]
            [clojure.core.reducers :as r])
  (:import [java.util Arrays]))

;; Enable optimizations for performance-critical code
(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn- local-maxima
  "Find local maxima in intensity data.
   Optimized version with option for primitive arrays for large datasets."
  [intensities window-size & {:keys [use-arrays?] :or {use-arrays? true}}]
  (if (and use-arrays? (> (count intensities) 1000))
    ;; Optimized path for large datasets
    (local-maxima-optimized intensities window-size)
    ;; Standard path for smaller datasets
    (local-maxima-standard intensities window-size)))

(defn- local-maxima-standard
  "Standard local maxima detection using collections."
  [intensities window-size]
  (let [n (count intensities)
        half-window (quot window-size 2)]
    (loop [i half-window
           peaks (transient [])]
      (if (>= i (- n half-window))
        (persistent! peaks)
        (let [center-val (nth intensities i)
              window-start (max 0 (- i half-window))
              window-end (min n (+ i half-window 1))
              is-maximum? (loop [j window-start
                                max? true]
                           (if (or (not max?) (>= j window-end))
                             max?
                             (recur (inc j) 
                                    (and max? (<= (nth intensities j) center-val)))))]
          (recur (inc i)
                 (if is-maximum?
                   (conj! peaks i)
                   peaks)))))))

(defn- local-maxima-optimized  
  "Optimized local maxima detection using primitive arrays."
  [intensities window-size]
  (let [intensity-array (if (instance? (Class/forName "[D") intensities)
                          intensities
                          (double-array intensities))
        n (alength ^doubles intensity-array)
        half-window (quot window-size 2)
        peaks (transient [])]
    (loop [i (long half-window)]
      (if (>= i (- n half-window))
        (persistent! peaks)
        (let [center-val (aget ^doubles intensity-array i)]
          (loop [j (long (max 0 (- i half-window)))
                 is-max true]
            (if (or (not is-max) (>= j (min n (+ i half-window 1))))
              (recur (inc i)
                     (if is-max
                       (conj! peaks i)
                       peaks))
              (recur (inc j)
                     (and is-max (<= (aget ^doubles intensity-array j) center-val)))))))))))

(defn- calculate-snr
  "Calculate signal-to-noise ratio for peak candidates.
   Optimized for batch processing."
  [intensities indices noise-estimate]
  (mapv (fn [idx]
          (let [signal (nth intensities idx)]
            (if (pos? noise-estimate)
              (/ signal noise-estimate)
              Double/MAX_VALUE)))
        indices))

(defn- estimate-noise
  "Estimate noise level using median absolute deviation.
   Optimized implementation with option for primitive arrays."
  [intensities & {:keys [use-arrays?] :or {use-arrays? true}}]
  (if (and use-arrays? (> (count intensities) 1000))
    (estimate-noise-optimized intensities)
    (estimate-noise-standard intensities)))

(defn- estimate-noise-standard
  "Standard noise estimation using collections."
  [intensities]
  (let [sorted (sort intensities)
        n (count sorted)
        median (nth sorted (quot n 2))
        deviations (persistent!
                    (reduce (fn [acc val]
                              (conj! acc (math/abs (- val median))))
                            (transient [])
                            intensities))
        mad (nth (sort deviations) (quot (count deviations) 2))]
    (* 1.4826 mad))) ; MAD to standard deviation conversion factor

(defn- estimate-noise-optimized
  "Optimized noise estimation using primitive arrays and sorting."
  [intensities]
  (let [intensity-array (if (instance? (Class/forName "[D") intensities)
                          intensities
                          (double-array intensities))
        n (alength ^doubles intensity-array)
        copy (double-array n)]
    (System/arraycopy intensity-array 0 copy 0 n)
    (Arrays/sort copy)
    (let [median (aget copy (quot n 2))
          deviations (double-array n)]
      (dotimes [i n]
        (aset deviations i (math/abs (- (aget ^doubles intensity-array i) median))))
      (Arrays/sort deviations)
      (* 1.4826 (aget deviations (quot n 2))))))

(defn detect-peaks
  "Detect peaks in MALDI spectrum with optimizations.
   
   Parameters:
   - mz-values: m/z values (vector)
   - intensities: intensity values (vector)
   - options: map with optional parameters
     - :snr-threshold (default 3.0)
     - :window-size (default 5)
     - :min-intensity (default 0.0)
     
   Returns: vector of maps with :mz, :intensity, :index, :snr"
  ([mz-values intensities]
   (detect-peaks mz-values intensities {}))
  ([mz-values intensities {:keys [snr-threshold window-size min-intensity use-optimizations?]
                            :or {snr-threshold 3.0
                                 window-size 5
                                 min-intensity 0.0
                                 use-optimizations? true}}]
   (let [n (count intensities)
         _ (when (not= n (count mz-values))
             (throw (ex-info "Mismatch between mz-values and intensities lengths"
                           {:mz-count (count mz-values)
                            :intensity-count n})))
         
         ;; Choose optimization level based on data size
         use-arrays? (and use-optimizations? (> n 1000))
         
         ;; Pre-filter by minimum intensity for efficiency
         valid-indices (if use-arrays?
                         ;; Optimized filtering for large datasets
                         (let [intensity-array (double-array intensities)]
                           (persistent!
                            (loop [i 0
                                   acc (transient [])]
                              (if (>= i n)
                                acc
                                (recur (inc i)
                                       (if (>= (aget intensity-array i) min-intensity)
                                         (conj! acc i)
                                         acc))))))
                         ;; Standard filtering for smaller datasets
                         (persistent!
                          (reduce-kv (fn [acc idx val]
                                       (if (>= val min-intensity)
                                         (conj! acc idx)
                                         acc))
                                     (transient [])
                                     (vec intensities))))
         
         ;; Find local maxima among valid points
         local-max-indices (local-maxima intensities window-size :use-arrays? use-arrays?)
         candidate-indices (filter (set valid-indices) local-max-indices)
         
         ;; Calculate noise estimate once
         noise-est (estimate-noise intensities :use-arrays? use-arrays?)
         
         ;; Calculate SNR for all candidates at once
         snr-values (calculate-snr intensities candidate-indices noise-est)
         
         ;; Filter by SNR threshold and build result
         peaks (persistent!
                (reduce (fn [acc [idx snr]]
                          (if (>= snr snr-threshold)
                            (conj! acc {:mz (nth mz-values idx)
                                       :intensity (nth intensities idx)
                                       :index idx
                                       :snr snr})
                            acc))
                        (transient [])
                        (map vector candidate-indices snr-values)))]
     
     ;; Sort by intensity (highest first)
     (sort-by :intensity #(compare %2 %1) peaks))))

(defn centroid-peaks
  "Apply centroiding to detected peaks for improved mass accuracy.
   Uses weighted average of peak region."
  [mz-values intensities peaks {:keys [centroid-window]
                                 :or {centroid-window 3}}]
  (let [half-window (quot centroid-window 2)
        n (count intensities)]
    (mapv (fn [{:keys [index] :as peak}]
            (let [start (max 0 (- index half-window))
                  end (min n (+ index half-window 1))
                  region-indices (range start end)
                  region-mz (mapv #(nth mz-values %) region-indices)
                  region-intensities (mapv #(nth intensities %) region-indices)
                  total-intensity (reduce + region-intensities)
                  weighted-mz (if (pos? total-intensity)
                                 (/ (reduce + (map * region-mz region-intensities))
                                    total-intensity)
                                 (:mz peak))]
              (assoc peak :mz weighted-mz :centroided true)))
          peaks)))

(defn filter-peaks
  "Filter peaks by various criteria with optimized batch processing."
  [peaks {:keys [min-snr max-peaks intensity-threshold]
          :or {min-snr 0.0
               max-peaks Integer/MAX_VALUE
               intensity-threshold 0.0}}]
  (->> peaks
       (filter #(>= (:snr %) min-snr))
       (filter #(>= (:intensity %) intensity-threshold))
       (sort-by :intensity #(compare %2 %1))
       (take max-peaks)))

;; Optimized batch processing functions

(defn detect-peaks-batch
  "Detect peaks across multiple spectra with optimized batch processing.
   Options:
   - :parallel? (default false) - use parallel processing for large batches
   - :chunk-size (default 100) - number of spectra per chunk for parallel processing
   Returns map of spectrum-id -> peaks"
  [spectra-map {:keys [parallel? chunk-size] 
                :or {parallel? false chunk-size 100}
                :as options}]
  
  (if (and parallel? (> (count spectra-map) chunk-size))
    ;; Parallel processing for large batches
    (->> spectra-map
         (r/map (fn [[spectrum-id {:keys [mz intensities]}]]
                  [spectrum-id (detect-peaks mz intensities options)]))
         (r/fold (partial merge-with concat) hash-map))
    
    ;; Sequential processing (optimized with transients)
    (persistent!
     (reduce-kv (fn [acc spectrum-id {:keys [mz intensities]}]
                  (assoc! acc spectrum-id 
                          (detect-peaks mz intensities options)))
                (transient {})
                spectra-map))))

(defn peak-statistics
  "Calculate statistics for peak detection results."
  [peaks]
  {:count (count peaks)
   :mean-intensity (if (seq peaks)
                     (/ (reduce + (map :intensity peaks)) (count peaks))
                     0.0)
   :median-intensity (if (seq peaks)
                       (nth (sort (map :intensity peaks))
                            (quot (count peaks) 2))
                       0.0)
   :mean-snr (if (seq peaks)
               (/ (reduce + (map :snr peaks)) (count peaks))
               0.0)
   :intensity-range [(apply min (map :intensity peaks))
                     (apply max (map :intensity peaks))]})

;; Performance monitoring and benchmarking

(defn benchmark-peak-detection
  "Benchmark peak detection with different optimization settings.
   Returns performance metrics."
  [mz-values intensities options]
  (let [base-options (dissoc options :use-optimizations?)
        
        ;; Test standard implementation
        start-time (System/nanoTime)
        standard-peaks (detect-peaks mz-values intensities 
                                    (assoc base-options :use-optimizations? false))
        standard-time (/ (- (System/nanoTime) start-time) 1000000.0)
        
        ;; Test optimized implementation
        start-time (System/nanoTime)
        optimized-peaks (detect-peaks mz-values intensities 
                                     (assoc base-options :use-optimizations? true))
        optimized-time (/ (- (System/nanoTime) start-time) 1000000.0)]
    
    {:standard {:time-ms standard-time
                :peak-count (count standard-peaks)}
     :optimized {:time-ms optimized-time
                 :peak-count (count optimized-peaks)}
     :speedup (if (pos? optimized-time) (/ standard-time optimized-time) 1.0)
     :data-size (count intensities)}))

(defn performance-profile
  "Create a performance profile for different data sizes."
  [& {:keys [sizes test-data-fn]
      :or {sizes [100 500 1000 5000 10000]
           test-data-fn (fn [size] 
                         [(range size) 
                          (repeatedly size #(+ 10 (rand 1000)))])}}]
  (mapv (fn [size]
          (let [[mz intensities] (test-data-fn size)
                benchmark (benchmark-peak-detection mz intensities {})]
            (assoc benchmark :size size)))
        sizes))

;; Memory usage utilities

(defn estimate-memory-usage
  "Estimate memory usage for peak detection on given data size."
  [data-size]
  (let [base-memory (* data-size 8 4) ; 4 double arrays
        temp-memory (* data-size 8 2) ; temporary arrays for sorting
        result-memory (* (/ data-size 100) 32)] ; estimated peak result size
    {:base-mb (/ base-memory 1024.0 1024.0)
     :temp-mb (/ temp-memory 1024.0 1024.0) 
     :result-mb (/ result-memory 1024.0 1024.0)
     :total-mb (/ (+ base-memory temp-memory result-memory) 1024.0 1024.0)}))

;; Export main functions
(def ^:export detect-peaks detect-peaks)
(def ^:export centroid-peaks centroid-peaks)
(def ^:export filter-peaks filter-peaks)
(def ^:export detect-peaks-batch detect-peaks-batch)
(def ^:export benchmark-peak-detection benchmark-peak-detection)
(def ^:export performance-profile performance-profile)(ns maldi-clj.peaks
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
