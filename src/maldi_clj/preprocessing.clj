(ns maldi-clj.preprocessing
  "High-performance preprocessing functions using dtype-next operations"
  (:require [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [tech.v3.dataset :as ds]))

;; ==============================================================================
;; HIGH-PERFORMANCE TRANSFORMATION METHODS - Using dtype-next operations
;; ==============================================================================

(defn sqrt-stabilize
  "Apply square root variance stabilization using efficient dtype-next operations"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use fast SIMD-optimized operations for dtype-next containers
    (dfn/sqrt (dfn/max intensities 0.0))
    ;; Fallback for legacy vectors
    (mapv #(Math/sqrt (max 0.0 %)) intensities)))

(defn log-transform
  "Apply natural logarithm transformation using efficient dtype-next operations"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use fast SIMD operations with proper handling of zero/negative values
    (dfn/log (dfn/max intensities Double/MIN_VALUE))
    ;; Fallback for legacy vectors
    (mapv #(Math/log (max % Double/MIN_VALUE)) intensities)))

(defn log10-transform
  "Apply base-10 logarithm transformation using efficient dtype-next operations"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use fast operations - dfn/log10 available in recent versions
    (dfn// (dfn/log (dfn/max intensities Double/MIN_VALUE)) (Math/log 10))
    ;; Fallback for legacy vectors
    (mapv #(Math/log10 (max % Double/MIN_VALUE)) intensities)))

(defn log2-transform
  "Apply base-2 logarithm transformation using efficient dtype-next operations"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use fast operations with proper base conversion
    (dfn// (dfn/log (dfn/max intensities Double/MIN_VALUE)) (Math/log 2))
    ;; Fallback for legacy vectors
    (mapv #(/ (Math/log (max % Double/MIN_VALUE)) (Math/log 2)) intensities)))

(defn transform-intensities
  "Apply various intensity transformations with automatic dtype-next optimization"
  [intensities method]
  (case method
    :sqrt (sqrt-stabilize intensities)
    :log (log-transform intensities)
    :log10 (log10-transform intensities)
    :log2 (log2-transform intensities)
    :none intensities
    (throw (ex-info "Unknown transformation method" {:method method}))))

;; ==============================================================================
;; HIGH-PERFORMANCE NORMALIZATION METHODS - Using dtype-next reductions
;; ==============================================================================

(defn tic-normalize
  "Normalize intensities by Total Ion Current using efficient operations"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use fast SIMD-optimized sum and division
    (let [tic (dfn/sum intensities)]
      (if (pos? tic)
        (dfn// intensities tic)
        intensities))
    ;; Fallback for legacy vectors
    (let [tic (reduce + intensities)]
      (if (pos? tic)
        (mapv #(/ % tic) intensities)
        intensities))))

(defn median-normalize
  "Normalize intensities by median value using efficient operations"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use dtype-next efficient median computation
    (let [median-val (if (zero? (dtype/ecount intensities))
                       1.0
                       (dfn/median intensities))]
      (if (zero? median-val)
        intensities
        (dfn// intensities median-val)))
    ;; Fallback for legacy vectors
    (let [median-val (if (empty? intensities)
                       1.0
                       (nth (sort intensities) (quot (count intensities) 2)))]
      (if (zero? median-val)
        intensities
        (mapv #(/ % median-val) intensities)))))

(defn calibrate-intensities
  "Apply various normalization methods with automatic optimization"
  [intensities method]
  (case method
    :tic (tic-normalize intensities)
    :median (median-normalize intensities)
    ;; PQN normalization would require multiple spectra for reference
    :none intensities
    (throw (ex-info "Unknown calibration method" {:method method}))))

;; ==============================================================================
;; HIGH-PERFORMANCE SMOOTHING FUNCTIONS - Using dtype-next windowing
;; ==============================================================================

(defn moving-average-smooth
  "Apply moving average smoothing using efficient dtype-next operations"
  [intensities half-window-size & {:keys [weighted?] :or {weighted? false}}]
  (if (dtype/reader? intensities)
    ;; Use dtype-next efficient operations
    (let [n (dtype/ecount intensities)
          window-size (inc (* 2 half-window-size))]
      (if weighted?
        ;; Weighted moving average using efficient operations
        (let [result (dtype/make-container :jvm-heap (dtype/elemwise-datatype intensities) n)]
          ;; Simple triangular weighting for now - could be optimized further
          (dotimes [i n]
            (let [start (max 0 (- i half-window-size))
                  end (min n (+ i half-window-size 1))
                  length (- end start)]
              (if (> length 0)
                (let [window (dtype/sub-buffer intensities start length)
                      center (- i start)
                      weights (dtype/make-container :jvm-heap :float64
                                                    (map #(- half-window-size (Math/abs (- % center)))
                                                         (range (dtype/ecount window))))
                      weighted-sum (dfn/sum (dfn/* window weights))
                      weight-sum (dfn/sum weights)]
                  (dtype/set-value! result i
                                    (if (zero? weight-sum) 0.0 (/ weighted-sum weight-sum))))
                (dtype/set-value! result i (dtype/get-value intensities i)))))
          result)
        ;; Simple moving average - implement manual rolling for compatibility
        (let [result (dtype/make-container :jvm-heap (dtype/elemwise-datatype intensities) n)]
          (dotimes [i n]
            (let [start (max 0 (- i half-window-size))
                  end (min n (+ i half-window-size 1))
                  length (- end start)]
              (if (> length 0)
                (let [window (dtype/sub-buffer intensities start length)]
                  (dtype/set-value! result i (dfn/mean window)))
                (dtype/set-value! result i (dtype/get-value intensities i)))))
          result)))
    ;; Fallback for legacy vectors
    (let [n (count intensities)]
      (mapv (fn [i]
              (let [start (max 0 (- i half-window-size))
                    end (min n (+ i half-window-size 1))
                    window (subvec intensities start end)]
                (if weighted?
                  ;; Weighted moving average (simple triangular weights)
                  (let [center (- i start)
                        weights (mapv #(- half-window-size (Math/abs (- % center))) (range (count window)))
                        weighted-sum (reduce + (map * window weights))
                        weight-sum (reduce + weights)]
                    (if (zero? weight-sum) 0.0 (/ weighted-sum weight-sum)))
                  ;; Simple moving average
                  (/ (reduce + window) (count window)))))
            (range n)))))

(defn smooth-intensities
  "Apply smoothing with automatic dtype-next optimization"
  [intensities method & {:keys [half-window-size] :or {half-window-size 10}}]
  (case method
    :moving-average (moving-average-smooth intensities half-window-size)
    :weighted-moving-average (moving-average-smooth intensities half-window-size :weighted? true)
    :none intensities
    (throw (ex-info "Unknown smoothing method" {:method method}))))

;; ==============================================================================
;; HIGH-PERFORMANCE SPECTRUM MANIPULATION - Using dtype-next filtering
;; ==============================================================================

(defn trim-spectrum
  "Trim spectrum to specified m/z range using efficient dtype-next operations"
  [spectrum min-mz max-mz]
  (let [mz-vals (:mz-values spectrum)
        intensities (:intensities spectrum)]

    (if (spectrum/dtype-spectrum? spectrum)
      ;; Use efficient dtype-next filtering
      (let [valid-mask (dfn/and (dfn/>= mz-vals min-mz)
                                (dfn/<= mz-vals max-mz))
            valid-indices (dtype/make-container :jvm-heap :int32
                                                (keep-indexed #(when %2 %1) (dtype/->vector valid-mask)))]
        (if (empty? valid-indices)
          (assoc spectrum
                 :mz-values (dtype/make-container :jvm-heap :float64 [])
                 :intensities (dtype/make-container :jvm-heap :float32 []))
          (let [trimmed-mz (dtype/indexed-buffer valid-indices mz-vals)
                trimmed-intensities (dtype/indexed-buffer valid-indices intensities)]
            (assoc spectrum
                   :mz-values trimmed-mz
                   :intensities trimmed-intensities))))

      ;; Fallback for legacy vectors
      (let [indices (keep-indexed
                     (fn [idx mz]
                       (when (<= min-mz mz max-mz) idx))
                     mz-vals)]
        (if (empty? indices)
          (assoc spectrum :mz-values [] :intensities [])
          (let [trimmed-mz (mapv #(nth mz-vals %) indices)
                trimmed-intensities (mapv #(nth intensities %) indices)]
            (assoc spectrum
                   :mz-values trimmed-mz
                   :intensities trimmed-intensities)))))))

(defn total-ion-current
  "Calculate Total Ion Current using efficient dtype-next operations"
  [spectrum]
  (let [intensities (:intensities spectrum)]
    (if (dtype/reader? intensities)
      (dfn/sum intensities)
      (reduce + 0 intensities))))

(defn is-empty-spectrum?
  "Check if spectrum is empty using efficient operations"
  [spectrum]
  (let [intensities (:intensities spectrum)]
    (if (dtype/reader? intensities)
      (zero? (dtype/ecount intensities))
      (empty? intensities))))

;; ==============================================================================
;; ENHANCED HIGH-PERFORMANCE PREPROCESSING PIPELINE
;; ==============================================================================

(defn preprocess-spectrum-enhanced
  "Enhanced preprocessing pipeline leveraging dtype-next for maximum performance"
  [spectrum & {:keys [transform-method smooth-method calibrate-method
                      smooth-window trim-range force-dtype-conversion?]
               :or {transform-method :none
                    smooth-method :none
                    calibrate-method :none
                    smooth-window 10
                    trim-range nil
                    force-dtype-conversion? true}}]
  (let [;; Convert to dtype-next format if needed and requested
        working-spectrum (if (and force-dtype-conversion? (spectrum/legacy-spectrum? spectrum))
                           (spectrum/clojure->spectrum spectrum)
                           spectrum)

        ;; Step 1: Trim if range specified  
        trimmed (if trim-range
                  (trim-spectrum working-spectrum (first trim-range) (second trim-range))
                  working-spectrum)

        ;; Step 2: Transform intensities using optimized operations
        transformed-intensities (transform-intensities (:intensities trimmed) transform-method)
        after-transform (assoc trimmed :intensities transformed-intensities)

        ;; Step 3: Smooth if requested using optimized operations
        smoothed-intensities (smooth-intensities
                              (:intensities after-transform)
                              smooth-method
                              :half-window-size smooth-window)
        after-smooth (assoc after-transform :intensities smoothed-intensities)

        ;; Step 4: Calibrate (normalize) using optimized operations
        calibrated-intensities (calibrate-intensities (:intensities after-smooth) calibrate-method)
        final-spectrum (assoc after-smooth :intensities calibrated-intensities)]

    ;; Validate result
    (if (spectrum/valid-spectrum? final-spectrum)
      final-spectrum
      (throw (ex-info "Enhanced preprocessing produced invalid spectrum"
                      {:original spectrum :result final-spectrum})))))

;; ==============================================================================
;; LEGACY COMPATIBILITY FUNCTIONS - For gradual migration
;; ==============================================================================

(defn preprocess-spectrum
  "Legacy preprocessing function with automatic dtype-next optimization
   Maintains API compatibility while providing performance benefits"
  [spectrum & {:keys [sqrt-stabilize? tic-normalize?]
               :or {sqrt-stabilize? false tic-normalize? false}}]
  (let [;; Automatically convert to dtype-next for processing
        efficient-spectrum (if (spectrum/legacy-spectrum? spectrum)
                             (spectrum/clojure->spectrum spectrum)
                             spectrum)

        ;; Apply legacy operations using new efficient methods
        intensities (:intensities efficient-spectrum)
        processed-intensities (cond-> intensities
                                sqrt-stabilize? sqrt-stabilize
                                tic-normalize? tic-normalize)
        new-spectrum (assoc efficient-spectrum :intensities processed-intensities)]

    ;; Validate the processed spectrum
    (if (spectrum/valid-spectrum? new-spectrum)
      ;; Return in same format as input for compatibility
      (if (spectrum/legacy-spectrum? spectrum)
        (spectrum/spectrum->clojure new-spectrum)
        new-spectrum)
      (throw (ex-info "Processed spectrum failed validation"
                      {:errors (spectrum/explain-spectrum-errors new-spectrum)
                       :original spectrum
                       :processed new-spectrum})))))

;; ==============================================================================
;; PERFORMANCE UTILITIES - For benchmarking and optimization
;; ==============================================================================

(defn benchmark-preprocessing
  "Benchmark preprocessing operations to demonstrate performance improvements"
  [n-points]
  (let [;; Create test data
        mz-data (range 100.0 (+ 100.0 n-points) 1.0)
        intensity-data (repeatedly n-points #(rand 10000))

        ;; Legacy spectrum
        legacy-spectrum (spectrum/create-spectrum-legacy "test" mz-data intensity-data {})

        ;; Efficient spectrum
        efficient-spectrum (spectrum/create-spectrum "test" mz-data intensity-data {})]

    (println (str "Benchmarking with " n-points " data points"))

    ;; Benchmark sqrt stabilization
    (println "\nSqrt stabilization:")
    (print "Legacy: ")
    (time (sqrt-stabilize (:intensities legacy-spectrum)))
    (print "Efficient: ")
    (time (sqrt-stabilize (:intensities efficient-spectrum)))

    ;; Benchmark TIC normalization
    (println "\nTIC normalization:")
    (print "Legacy: ")
    (time (tic-normalize (:intensities legacy-spectrum)))
    (print "Efficient: ")
    (time (tic-normalize (:intensities efficient-spectrum)))

    ;; Benchmark full preprocessing
    (println "\nFull preprocessing pipeline:")
    (print "Legacy: ")
    (time (preprocess-spectrum legacy-spectrum :sqrt-stabilize? true :tic-normalize? true))
    (print "Efficient: ")
    (time (preprocess-spectrum-enhanced efficient-spectrum
                                        :transform-method :sqrt
                                        :calibrate-method :tic))

    nil))