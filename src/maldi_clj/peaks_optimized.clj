(ns maldi-clj.peaks-optimized
  "Highly optimized peak detection algorithms for MALDI mass spectrometry data.
   Advanced optimizations including parallel processing, memory-efficient algorithms,
   and specialized data structures for large-scale analysis."
  (:require [clojure.math :as math]
            [clojure.core.reducers :as r])
  (:import [java.util.concurrent ForkJoinPool]
           [java.util Arrays]))

;; Type hints and optimized data structures

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn- array-local-maxima
  "Ultra-fast local maxima detection using primitive arrays."
  [^doubles intensities ^long window-size]
  (let [n (alength intensities)
        half-window (quot window-size 2)
        peaks (transient [])]
    (loop [i (long half-window)]
      (if (>= i (- n half-window))
        (persistent! peaks)
        (let [center-val (aget intensities i)]
          (loop [j (long (max 0 (- i half-window)))
                 is-max true]
            (if (or (not is-max) (>= j (min n (+ i half-window 1))))
              (recur (inc i)
                     (if is-max
                       (conj! peaks i)
                       peaks))
              (recur (inc j)
                     (and is-max (<= (aget intensities j) center-val))))))))))

(defn- fast-noise-estimate
  "Optimized noise estimation using primitive arrays and efficient sorting."
  [^doubles intensities]
  (let [n (alength intensities)
        copy (double-array n)]
    (System/arraycopy intensities 0 copy 0 n)
    (Arrays/sort copy)
    (let [median (aget copy (quot n 2))
          deviations (double-array n)]
      (dotimes [i n]
        (aset deviations i (math/abs (- (aget intensities i) median))))
      (Arrays/sort deviations)
      (* 1.4826 (aget deviations (quot n 2))))))

(defn- parallel-peak-detection
  "Parallel peak detection for large datasets using core.reducers."
  [mz-values intensities options]
  (let [{:keys [snr-threshold window-size min-intensity chunk-size]
         :or {snr-threshold 3.0
              window-size 5
              min-intensity 0.0
              chunk-size 10000}} options
        n (count intensities)
        intensity-array (double-array intensities)
        mz-array (double-array mz-values)]
    
    ;; Process in chunks for better memory usage and parallelization
    (if (< n chunk-size)
      ;; Small dataset - use regular algorithm
      (array-local-maxima intensity-array window-size)
      
      ;; Large dataset - use parallel processing
      (let [chunks (partition-all chunk-size (range n))
            noise-est (fast-noise-estimate intensity-array)]
        
        ;; Process chunks in parallel
        (->> chunks
             (r/map (fn [chunk-indices]
                      (let [chunk-start (first chunk-indices)
                            chunk-end (inc (last chunk-indices))
                            chunk-size (- chunk-end chunk-start)
                            chunk-intensities (double-array chunk-size)
                            chunk-mz (double-array chunk-size)]
                        
                        ;; Copy chunk data
                        (dotimes [i chunk-size]
                          (let [global-idx (+ chunk-start i)]
                            (aset chunk-intensities i (aget intensity-array global-idx))
                            (aset chunk-mz i (aget mz-array global-idx))))
                        
                        ;; Find peaks in chunk
                        (let [local-peaks (array-local-maxima chunk-intensities window-size)]
                          ;; Convert to global indices and filter by criteria
                          (->> local-peaks
                               (map #(+ % chunk-start))
                               (filter (fn [idx]
                                         (let [intensity (aget intensity-array idx)
                                               snr (if (pos? noise-est)
                                                     (/ intensity noise-est)
                                                     Double/MAX_VALUE)]
                                           (and (>= intensity min-intensity)
                                                (>= snr snr-threshold)))))
                               (map (fn [idx]
                                      {:mz (aget mz-array idx)
                                       :intensity (aget intensity-array idx)
                                       :index idx
                                       :snr (/ (aget intensity-array idx) noise-est)})))))))
             (r/fold concat [])
             (sort-by :intensity #(compare %2 %1)))))))

(defn- memory-efficient-centroiding
  "Memory-efficient centroiding using streaming computation."
  [^doubles mz-array ^doubles intensity-array peaks centroid-window]
  (let [half-window (quot centroid-window 2)
        n (alength intensity-array)]
    (mapv (fn [{:keys [index] :as peak}]
            (let [start (max 0 (- ^long index half-window))
                  end (min n (+ ^long index half-window 1))
                  
                  ;; Stream computation to avoid intermediate collections
                  [total-weight weighted-sum] 
                  (loop [i start
                         total 0.0
                         weighted 0.0]
                    (if (>= i end)
                      [total weighted]
                      (let [intensity (aget intensity-array i)
                            mz (aget mz-array i)]
                        (recur (inc i)
                               (+ total intensity)
                               (+ weighted (* intensity mz))))))
                  
                  centroided-mz (if (pos? total-weight)
                                  (/ weighted-sum total-weight)
                                  (:mz peak))]
              (assoc peak :mz centroided-mz :centroided true)))
          peaks)))

(defn detect-peaks-optimized
  "Highly optimized peak detection with advanced performance features.
   
   Features:
   - Primitive array operations for maximum speed
   - Parallel processing for large datasets  
   - Memory-efficient algorithms
   - Configurable chunking for optimal memory usage
   
   Additional options:
   - :parallel? (default true) - enable parallel processing
   - :chunk-size (default 10000) - chunk size for parallel processing
   - :use-primitive-arrays? (default true) - use primitive arrays"
  ([mz-values intensities]
   (detect-peaks-optimized mz-values intensities {}))
  ([mz-values intensities {:keys [parallel? use-primitive-arrays?]
                            :or {parallel? true
                                 use-primitive-arrays? true}
                            :as options}]
   
   (cond
     ;; Ultra-fast path for primitive arrays and parallel processing
     (and parallel? use-primitive-arrays? (> (count intensities) 1000))
     (parallel-peak-detection mz-values intensities options)
     
     ;; Fast path for primitive arrays, single-threaded
     use-primitive-arrays?
     (let [intensity-array (double-array intensities)
           mz-array (double-array mz-values)
           {:keys [snr-threshold window-size min-intensity]
            :or {snr-threshold 3.0 window-size 5 min-intensity 0.0}} options
           
           peak-indices (array-local-maxima intensity-array window-size)
           noise-est (fast-noise-estimate intensity-array)]
       
       (->> peak-indices
            (filter (fn [idx]
                      (let [intensity (aget intensity-array idx)
                            snr (if (pos? noise-est)
                                  (/ intensity noise-est)  
                                  Double/MAX_VALUE)]
                        (and (>= intensity min-intensity)
                             (>= snr snr-threshold)))))
            (map (fn [idx]
                   {:mz (aget mz-array idx)
                    :intensity (aget intensity-array idx)
                    :index idx
                    :snr (/ (aget intensity-array idx) noise-est)}))
            (sort-by :intensity #(compare %2 %1))))
     
     ;; Fallback to standard implementation
     :else
     (throw (ex-info "Standard peak detection not implemented in optimized module"
                   {:suggestion "Use maldi-clj.peaks/detect-peaks for standard implementation"})))))

(defn centroid-peaks-optimized
  "Optimized centroiding using primitive arrays and streaming computation."
  [mz-values intensities peaks options]
  (let [{:keys [centroid-window use-primitive-arrays?]
         :or {centroid-window 3 use-primitive-arrays? true}} options]
    
    (if use-primitive-arrays?
      (memory-efficient-centroiding (double-array mz-values)
                                    (double-array intensities)
                                    peaks
                                    centroid-window)
      (throw (ex-info "Non-primitive array centroiding not implemented in optimized module"
                    {:suggestion "Use maldi-clj.peaks/centroid-peaks for standard implementation"})))))

(defn benchmark-peak-detection
  "Benchmark different peak detection implementations."
  [mz-values intensities options]
  (let [standard-time (time (doall (require 'maldi-clj.peaks)
                                   ((resolve 'maldi-clj.peaks/detect-peaks) 
                                    mz-values intensities options)))
        
        optimized-time (time (doall (detect-peaks-optimized mz-values intensities options)))
        
        parallel-time (time (doall (detect-peaks-optimized mz-values intensities 
                                                          (assoc options :parallel? true))))]
    
    {:standard-ms standard-time
     :optimized-ms optimized-time  
     :parallel-ms parallel-time
     :speedup-optimized (/ standard-time optimized-time)
     :speedup-parallel (/ standard-time parallel-time)}))

;; Batch processing with advanced optimizations

(defn detect-peaks-batch-optimized
  "Batch peak detection with memory pooling and work stealing."
  [spectra-map options]
  (let [{:keys [max-threads]
         :or {max-threads (.availableProcessors (Runtime/getRuntime))}} options
        
        pool (ForkJoinPool. max-threads)]
    
    (try
      (->> spectra-map
           (r/map (fn [[spectrum-id {:keys [mz intensities]}]]
                    [spectrum-id (detect-peaks-optimized mz intensities options)]))
           (r/fold (partial merge-with concat) hash-map))
      
      (finally
        (.shutdown pool)))))

;; Export optimized functions
(def ^:export detect-peaks-optimized detect-peaks-optimized)
(def ^:export centroid-peaks-optimized centroid-peaks-optimized)
(def ^:export benchmark-peak-detection benchmark-peak-detection)(ns maldi-clj.peaks-optimized
  "Ultra-high-performance peak detection using dtype-next buffers"
  (:require [maldi-clj.spectrum :as spectrum]
            [maldi-clj.baseline :as baseline]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]))

;; =============================================================================
;; High-Performance Peak Data Structures
;; =============================================================================

(defn create-peak-buffers
  "Create dtype-next buffers for peak data instead of maps/vectors"
  [capacity]
  {:mz-buffer (dtype/make-container :jvm-heap :float64 capacity)
   :intensity-buffer (dtype/make-container :jvm-heap :float32 capacity)
   :snr-buffer (dtype/make-container :jvm-heap :float32 capacity)
   :noise-buffer (dtype/make-container :jvm-heap :float32 capacity)
   :index-buffer (dtype/make-container :jvm-heap :int32 capacity)
   :count 0})

(defn add-peak-to-buffers!
  "Add a peak to the buffer structure - much faster than conj"
  [peak-buffers idx mz intensity snr noise]
  (let [current-count (:count peak-buffers)]
    (dtype/set-value! (:mz-buffer peak-buffers) current-count mz)
    (dtype/set-value! (:intensity-buffer peak-buffers) current-count intensity)
    (dtype/set-value! (:snr-buffer peak-buffers) current-count snr)
    (dtype/set-value! (:noise-buffer peak-buffers) current-count noise)
    (dtype/set-value! (:index-buffer peak-buffers) current-count idx)
    (assoc peak-buffers :count (inc current-count))))

(defn truncate-peak-buffers
  "Truncate buffers to actual size - removes unused capacity"
  [peak-buffers]
  (let [actual-count (:count peak-buffers)]
    (if (zero? actual-count)
      peak-buffers
      {:mz-buffer (dtype/sub-buffer (:mz-buffer peak-buffers) 0 actual-count)
       :intensity-buffer (dtype/sub-buffer (:intensity-buffer peak-buffers) 0 actual-count)
       :snr-buffer (dtype/sub-buffer (:snr-buffer peak-buffers) 0 actual-count)
       :noise-buffer (dtype/sub-buffer (:noise-buffer peak-buffers) 0 actual-count)
       :index-buffer (dtype/sub-buffer (:index-buffer peak-buffers) 0 actual-count)
       :count actual-count})))

;; =============================================================================
;; Optimized Local Maxima Detection
;; =============================================================================

(defn find-local-maxima-optimized
  "Ultra-fast local maxima detection using dtype-next operations"
  [intensities & {:keys [window-half-size min-intensity]
                  :or {window-half-size 1 min-intensity 0}}]
  (if (dtype/reader? intensities)
    (let [n (dtype/ecount intensities)
          ;; Pre-allocate result buffer - much faster than building vector
          maxima-buffer (dtype/make-container :jvm-heap :int32 (quot n 10)) ; Estimate ~10% are maxima
          maxima-count (atom 0)]

      (dotimes [i n]
        (let [current-intensity (dtype/get-value intensities i)]
          (when (>= current-intensity min-intensity)
            (let [start (max 0 (- i window-half-size))
                  end (min n (+ i window-half-size 1))
                  ;; Use vectorized max instead of loop
                  window (dtype/sub-buffer intensities start (- end start))
                  max-in-window (dfn/reduce-max window)]

              (when (= current-intensity max-in-window)
                (dtype/set-value! maxima-buffer @maxima-count i)
                (swap! maxima-count inc))))))

      ;; Return properly sized buffer
      (dtype/sub-buffer maxima-buffer 0 @maxima-count))

    ;; Fallback for non-dtype data
    (find-local-maxima intensities :window-half-size window-half-size :min-intensity min-intensity)))

;; =============================================================================
;; Vectorized SNR Calculation
;; =============================================================================

(defn calculate-snr-vectorized
  "Vectorized SNR calculation for all peaks at once"
  [intensities noise-estimates indices]
  (if (dtype/reader? intensities)
    (let [n-peaks (dtype/ecount indices)
          peak-intensities (dtype/make-container :jvm-heap :float32 n-peaks)
          peak-noises (dtype/make-container :jvm-heap :float32 n-peaks)]

      ;; Extract peak intensities and noises vectorized
      (dotimes [i n-peaks]
        (let [idx (dtype/get-value indices i)]
          (dtype/set-value! peak-intensities i (dtype/get-value intensities idx))
          (dtype/set-value! peak-noises i
                            (if (number? noise-estimates)
                              noise-estimates
                              (dtype/get-value noise-estimates idx)))))

      ;; Vectorized SNR calculation: (intensity - noise) / noise
      (dfn/max (dfn// (dfn/- peak-intensities peak-noises) peak-noises) 0.0))

    ;; Fallback
    (mapv #(let [intensity (nth intensities %)
                 noise (if (number? noise-estimates) noise-estimates (nth noise-estimates %))]
             (max (/ (- intensity noise) noise) 0.0))
          indices)))

;; =============================================================================
;; Ultra-Fast Peak Detection
;; =============================================================================

(defn detect-peaks-optimized
  "Ultra-high-performance peak detection using dtype-next buffers throughout"
  [spectrum & {:keys [snr-threshold window-half-size noise-method min-intensity baseline-corrected?]
               :or {snr-threshold 3.0 window-half-size 10 noise-method :mad
                    min-intensity 0 baseline-corrected? false}}]

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

        ;; Step 2: Optimized noise estimation
        noise-estimates (case noise-method
                          :mad (estimate-noise-mad intensities)
                          :local-mad (estimate-noise-local-mad intensities :window-half-size 25)
                          (throw (ex-info "Unknown noise method" {:method noise-method})))

        ;; Step 3: Ultra-fast local maxima detection
        maxima-indices (find-local-maxima-optimized intensities
                                                    :window-half-size window-half-size
                                                    :min-intensity min-intensity)

        ;; Step 4: Vectorized SNR calculation for all peaks at once
        snr-values (calculate-snr-vectorized intensities noise-estimates maxima-indices)

        ;; Step 5: Vectorized filtering by SNR threshold
        snr-mask (dfn/>= snr-values snr-threshold)
        n-maxima (dtype/ecount maxima-indices)

        ;; Count peaks that pass SNR filter
        peak-count (dtype/ecount (dtype/select snr-mask [true]))

        ;; Step 6: Build optimized peak buffers
        peak-buffers (create-peak-buffers peak-count)

        ;; Fill buffers efficiently
        final-buffers
        (loop [i 0
               buffers peak-buffers]
          (if (>= i n-maxima)
            buffers
            (if (dtype/get-value snr-mask i)
              (let [idx (dtype/get-value maxima-indices i)
                    mz (dtype/get-value mz-values idx)
                    intensity (dtype/get-value intensities idx)
                    snr (dtype/get-value snr-values i)
                    noise (if (number? noise-estimates)
                            noise-estimates
                            (dtype/get-value noise-estimates idx))]
                (recur (inc i) (add-peak-to-buffers! buffers idx mz intensity snr noise)))
              (recur (inc i) buffers))))]

    ;; Return optimized peak structure
    {:id (:id spectrum)
     :peak-buffers (truncate-peak-buffers final-buffers)
     :metadata (merge (:metadata spectrum)
                      {:peak-detection {:snr-threshold snr-threshold
                                        :window-half-size window-half-size
                                        :noise-method noise-method
                                        :baseline-corrected? baseline-corrected?
                                        :total-peaks (:count final-buffers)
                                        :total-maxima n-maxima
                                        :optimized-buffers? true}})}))

;; =============================================================================
;; Buffer-based Peak Analysis (replaces mapv operations)
;; =============================================================================

(defn get-peak-intensities-fast
  "Ultra-fast peak intensity extraction - no mapv!"
  [optimized-peaks]
  (get-in optimized-peaks [:peak-buffers :intensity-buffer]))

(defn get-peak-mz-values-fast
  "Ultra-fast peak m/z extraction - no mapv!"
  [optimized-peaks]
  (get-in optimized-peaks [:peak-buffers :mz-buffer]))

(defn get-peak-snr-values-fast
  "Ultra-fast peak SNR extraction - no mapv!"
  [optimized-peaks]
  (get-in optimized-peaks [:peak-buffers :snr-buffer]))

(defn filter-peaks-by-snr-vectorized
  "Ultra-fast vectorized peak filtering"
  [optimized-peaks min-snr]
  (let [buffers (:peak-buffers optimized-peaks)
        snr-buffer (:snr-buffer buffers)
        mask (dfn/>= snr-buffer min-snr)
        filtered-indices (dtype/select mask [true])
        n-filtered (dtype/ecount filtered-indices)]

    (if (zero? n-filtered)
      (assoc optimized-peaks :peak-buffers (create-peak-buffers 0))
      (let [new-buffers (create-peak-buffers n-filtered)]
        ;; Copy filtered peaks using vectorized operations
        (dotimes [i n-filtered]
          (let [original-idx (dtype/get-value filtered-indices i)]
            (add-peak-to-buffers! new-buffers
                                  (dtype/get-value (:index-buffer buffers) original-idx)
                                  (dtype/get-value (:mz-buffer buffers) original-idx)
                                  (dtype/get-value (:intensity-buffer buffers) original-idx)
                                  (dtype/get-value (:snr-buffer buffers) original-idx)
                                  (dtype/get-value (:noise-buffer buffers) original-idx))))
        (assoc optimized-peaks :peak-buffers (truncate-peak-buffers new-buffers))))))

(defn peak-statistics-vectorized
  "Ultra-fast vectorized peak statistics"
  [optimized-peaks]
  (let [buffers (:peak-buffers optimized-peaks)
        intensities (:intensity-buffer buffers)
        snr-values (:snr-buffer buffers)
        n-peaks (:count buffers)]

    (when (pos? n-peaks)
      {:count n-peaks
       :intensity-stats {:min (dfn/reduce-min intensities)
                         :max (dfn/reduce-max intensities)
                         :mean (dfn/mean intensities)
                         :sum (dfn/sum intensities)}
       :snr-stats {:min (dfn/reduce-min snr-values)
                   :max (dfn/reduce-max snr-values)
                   :mean (dfn/mean snr-values)
                   :sum (dfn/sum snr-values)}})))

;; =============================================================================
;; Compatibility Layer (for existing API)  
;; =============================================================================

(defn optimized-peaks->legacy
  "Convert optimized peak buffers back to legacy map format for compatibility"
  [optimized-peaks]
  (let [buffers (:peak-buffers optimized-peaks)
        n-peaks (:count buffers)]
    {:id (:id optimized-peaks)
     :peaks (mapv (fn [i]
                    {:mz (dtype/get-value (:mz-buffer buffers) i)
                     :intensity (dtype/get-value (:intensity-buffer buffers) i)
                     :snr (dtype/get-value (:snr-buffer buffers) i)
                     :noise (dtype/get-value (:noise-buffer buffers) i)
                     :index (dtype/get-value (:index-buffer buffers) i)})
                  (range n-peaks))
     :metadata (:metadata optimized-peaks)}))