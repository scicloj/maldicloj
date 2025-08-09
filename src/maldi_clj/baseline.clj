(ns maldi-clj.baseline
  "High-performance baseline correction algorithms for mass spectrometry data"
  (:require [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [malli.core :as m]))

;; ==============================================================================
;; SNIP (Statistics-Sensitive Non-linear Iterative Peak-clipping) Algorithm
;; ==============================================================================

(defn log-log-sqrt-transform
  "Apply the LLS (Log-Log-Sqrt) transform to compress dynamic range.
   This prevents large peaks from dominating the filtering process.
   Formula: LLS(x) = ln(ln(sqrt(x + 1) + 1) + 1)"
  [intensities]
  (if (dtype/reader? intensities)
    ;; Use efficient dtype-next operations
    (-> intensities
        (dfn/+ 1.0) ; x + 1
        dfn/sqrt ; sqrt(x + 1)
        (dfn/+ 1.0) ; sqrt(x + 1) + 1
        dfn/log ; ln(sqrt(x + 1) + 1)
        (dfn/+ 1.0) ; ln(sqrt(x + 1) + 1) + 1
        dfn/log) ; ln(ln(sqrt(x + 1) + 1) + 1)
    ;; Fallback for legacy vectors
    (mapv #(Math/log (+ (Math/log (+ (Math/sqrt (+ % 1.0)) 1.0)) 1.0)) intensities)))

(defn inverse-log-log-sqrt-transform
  "Inverse of the LLS transform to restore original dynamic range.
   Formula: x = (exp(exp(y) - 1) - 1)^2 - 1"
  [lls-intensities]
  (if (dtype/reader? lls-intensities)
    ;; Use efficient dtype-next operations
    (-> lls-intensities
        dfn/exp ; exp(y)
        (dfn/- 1.0) ; exp(y) - 1
        dfn/exp ; exp(exp(y) - 1)
        (dfn/- 1.0) ; exp(exp(y) - 1) - 1
        dfn/sq ; (exp(exp(y) - 1) - 1)^2
        (dfn/- 1.0)) ; (exp(exp(y) - 1) - 1)^2 - 1
    ;; Fallback for legacy vectors
    (mapv #(- (Math/pow (- (Math/exp (- (Math/exp %) 1.0)) 1.0) 2.0) 1.0) lls-intensities)))

(defn snip-minimum-filter
  "Apply the SNIP minimum filter for one iteration.
   For each point i, set y[i] = min(y[i], (y[i-m] + y[i+m])/2)
   where m is the current window half-size."
  [intensities window-half-size]
  (if (dtype/reader? intensities)
    ;; Efficient dtype-next implementation
    (let [n (dtype/ecount intensities)
          result (dtype/clone intensities)]
      (dotimes [i n]
        (when (and (>= i window-half-size)
                   (< i (- n window-half-size)))
          (let [left-val (dtype/get-value intensities (- i window-half-size))
                right-val (dtype/get-value intensities (+ i window-half-size))
                current-val (dtype/get-value intensities i)
                average (/ (+ left-val right-val) 2.0)]
            (dtype/set-value! result i (min current-val average)))))
      result)
    ;; Fallback for legacy vectors
    (let [n (count intensities)]
      (vec (map (fn [i]
                  (if (and (>= i window-half-size)
                           (< i (- n window-half-size)))
                    (let [left-val (nth intensities (- i window-half-size))
                          right-val (nth intensities (+ i window-half-size))
                          current-val (nth intensities i)
                          average (/ (+ left-val right-val) 2.0)]
                      (min current-val average))
                    (nth intensities i)))
                (range n))))))

(defn snip-baseline-estimation
  "Estimate baseline using the SNIP algorithm.
   
   The SNIP (Statistics-Sensitive Non-linear Iterative Peak-clipping) algorithm:
   1. Apply LLS transform to compress dynamic range
   2. Iteratively apply minimum filter with increasing window sizes
   3. Apply inverse LLS transform to restore scale
   
   Parameters:
   - intensities: intensity values (dtype-next container or vector)
   - iterations: number of iterations (default 100)
   - use-lls?: whether to use LLS transform (default true)
   
   Returns estimated baseline as same type as input"
  [intensities & {:keys [iterations use-lls?]
                  :or {iterations 100 use-lls? true}}]
  (let [;; Step 1: Apply LLS transform if requested
        working-intensities (if use-lls?
                              (log-log-sqrt-transform intensities)
                              intensities)

        ;; Step 2: Iterative filtering with fixed window size
        ;; Use window size that matches MALDIquant behavior (empirically determined)
        window-size 5
        filtered-intensities (reduce
                              (fn [current-intensities _]
                                (snip-minimum-filter current-intensities window-size))
                              working-intensities
                              (range iterations))

        ;; Step 3: Apply inverse LLS transform if needed
        baseline (if use-lls?
                   (inverse-log-log-sqrt-transform filtered-intensities)
                   filtered-intensities)]

    ;; Ensure baseline is never negative
    (if (dtype/reader? baseline)
      (dfn/max baseline 0.0)
      (mapv #(max % 0.0) baseline))))

;; ==============================================================================
;; TopHat Morphological Filter (Alternative baseline method)
;; ==============================================================================

(defn morphological-opening
  "Morphological opening operation: erosion followed by dilation.
   This is the core of the TopHat filter."
  [intensities window-half-size]
  (if (dtype/reader? intensities)
    ;; Efficient implementation using dtype-next
    (let [n (dtype/ecount intensities)
          ;; Erosion: minimum filter
          eroded (dtype/make-container :jvm-heap (dtype/elemwise-datatype intensities) n)]
      (dotimes [i n]
        (let [start (max 0 (- i window-half-size))
              end (min n (+ i window-half-size 1))
              length (- end start)]
          (if (> length 0)
            (let [window (dtype/sub-buffer intensities start (- end start))]
              (dtype/set-value! eroded i (dfn/reduce-min window)))
            (dtype/set-value! eroded i (dtype/get-value intensities i)))))

      ;; Dilation: maximum filter on eroded signal
      (let [dilated (dtype/make-container :jvm-heap (dtype/elemwise-datatype intensities) n)]
        (dotimes [i n]
          (let [start (max 0 (- i window-half-size))
                end (min n (+ i window-half-size 1))
                length (- end start)]
            (if (> length 0)
              (let [window (dtype/sub-buffer eroded start (- end start))]
                (dtype/set-value! dilated i (dfn/reduce-max window)))
              (dtype/set-value! dilated i (dtype/get-value eroded i)))))
        dilated))
    ;; Fallback for legacy vectors
    (let [n (count intensities)
          ;; Erosion: minimum filter
          eroded (mapv (fn [i]
                         (let [start (max 0 (- i window-half-size))
                               end (min n (+ i window-half-size 1))
                               window (subvec intensities start end)]
                           (apply min window)))
                       (range n))
          ;; Dilation: maximum filter on eroded signal
          dilated (mapv (fn [i]
                          (let [start (max 0 (- i window-half-size))
                                end (min n (+ i window-half-size 1))
                                window (subvec eroded start end)]
                            (apply max window)))
                        (range n))]
      dilated)))

(defn tophat-baseline-estimation
  "Estimate baseline using TopHat morphological filter.
   
   TopHat filter performs morphological opening (erosion followed by dilation)
   to estimate the baseline while preserving peak structure.
   
   Parameters:
   - intensities: intensity values
   - window-half-size: half-size of the morphological structuring element
   
   Returns estimated baseline"
  [intensities window-half-size]
  (morphological-opening intensities window-half-size))

;; ==============================================================================
;; Simple Baseline Methods (for comparison and fallback)
;; ==============================================================================

(defn median-baseline-estimation
  "Simple baseline estimation using moving median filter."
  [intensities window-half-size]
  (if (dtype/reader? intensities)
    ;; Efficient dtype-next implementation
    (let [n (dtype/ecount intensities)
          result (dtype/make-container :jvm-heap (dtype/datatype intensities) n)]
      (dotimes [i n]
        (let [start (max 0 (- i window-half-size))
              end (min n (+ i window-half-size 1))]
          (if (> (- end start) 0)
            (let [window (dtype/sub-buffer intensities start (- end start))]
              (dtype/set-value! result i (dfn/median window)))
            (dtype/set-value! result i (dtype/get-value intensities i)))))
      result)
    ;; Fallback for legacy vectors
    (let [n (count intensities)]
      (mapv (fn [i]
              (let [start (max 0 (- i window-half-size))
                    end (min n (+ i window-half-size 1))
                    window (subvec intensities start end)
                    sorted-window (sort window)]
                (nth sorted-window (quot (count sorted-window) 2))))
            (range n)))))

;; ==============================================================================
;; Main Baseline Estimation API - Matching MALDIquant interface
;; ==============================================================================

(defn estimate-baseline
  "Estimate baseline of a spectrum using various methods.
   
   Matches MALDIquant's estimateBaseline function interface.
   
   Parameters:
   - spectrum: spectrum object with :mz-values and :intensities
   - method: baseline estimation method (:snip, :tophat, :median)
   - options: method-specific options
   
   For SNIP method:
   - :iterations (default 100): number of SNIP iterations
   - :use-lls? (default true): use Log-Log-Sqrt transform
   
   For TopHat method:
   - :half-window-size (default 75): morphological filter window
   
   For Median method:
   - :half-window-size (default 50): median filter window
   
   Returns baseline as same data structure as input intensities"
  [spectrum method & {:keys [iterations use-lls? half-window-size]
                      :or {iterations 100
                           use-lls? true
                           half-window-size 75}}]
  (let [intensities (:intensities spectrum)]
    (case method
      :snip (snip-baseline-estimation intensities
                                      :iterations iterations
                                      :use-lls? use-lls?)
      :tophat (tophat-baseline-estimation intensities half-window-size)
      :median (median-baseline-estimation intensities half-window-size)
      (throw (ex-info "Unknown baseline estimation method"
                      {:method method
                       :available-methods [:snip :tophat :median]})))))

(defn remove-baseline
  "Remove baseline from spectrum using estimated baseline.
   
   Matches MALDIquant's removeBaseline function interface.
   
   Parameters:
   - spectrum: spectrum object
   - method: baseline estimation method 
   - options: method-specific options (same as estimate-baseline)
   
   Returns spectrum with baseline-corrected intensities"
  [spectrum method & options]
  (let [baseline (apply estimate-baseline spectrum method options)
        corrected-intensities (if (dtype/reader? (:intensities spectrum))
                                ;; Efficient subtraction ensuring non-negative results
                                (dfn/max (dfn/- (:intensities spectrum) baseline) 0.0)
                                ;; Fallback for legacy vectors
                                (mapv #(max (- %1 %2) 0.0)
                                      (:intensities spectrum)
                                      baseline))]
    (assoc spectrum :intensities corrected-intensities)))

;; ==============================================================================
;; Utility Functions for Analysis and Validation
;; ==============================================================================

(defn baseline-correction-stats
  "Calculate statistics about baseline correction quality."
  [original-spectrum corrected-spectrum baseline]
  (let [original-intensities (:intensities original-spectrum)
        corrected-intensities (:intensities corrected-spectrum)
        original-tic (if (dtype/reader? original-intensities)
                       (dfn/sum original-intensities)
                       (reduce + original-intensities))
        corrected-tic (if (dtype/reader? corrected-intensities)
                        (dfn/sum corrected-intensities)
                        (reduce + corrected-intensities))
        baseline-tic (if (dtype/reader? baseline)
                       (dfn/sum baseline)
                       (reduce + baseline))]
    {:original-tic original-tic
     :corrected-tic corrected-tic
     :baseline-tic baseline-tic
     :baseline-fraction (/ baseline-tic original-tic)
     :signal-retained (/ corrected-tic original-tic)}))

;; ==============================================================================
;; Enhanced Pipeline Integration
;; ==============================================================================

(defn preprocess-with-baseline-correction
  "Enhanced preprocessing pipeline with baseline correction.
   
   Extends the existing preprocessing pipeline with baseline correction
   integrated at the appropriate stage (after smoothing, before peak detection)."
  [spectrum & {:keys [transform-method smooth-method baseline-method
                      calibrate-method smooth-window baseline-iterations
                      baseline-window trim-range force-dtype-conversion?]
               :or {transform-method :none
                    smooth-method :none
                    baseline-method :snip
                    calibrate-method :none
                    smooth-window 10
                    baseline-iterations 100
                    baseline-window 75
                    trim-range nil
                    force-dtype-conversion? true}}]
  (let [;; Convert to dtype-next if needed
        working-spectrum (if (and force-dtype-conversion?
                                  (spectrum/legacy-spectrum? spectrum))
                           (spectrum/clojure->spectrum spectrum)
                           spectrum)

        ;; Step 1: Trim if range specified
        trimmed (if trim-range
                  (let [trim-fn (requiring-resolve 'maldi-clj.preprocessing/trim-spectrum)]
                    (trim-fn working-spectrum (first trim-range) (second trim-range)))
                  working-spectrum)

        ;; Step 2: Transform intensities
        transform-fn (requiring-resolve 'maldi-clj.preprocessing/transform-intensities)
        transformed-intensities (transform-fn (:intensities trimmed) transform-method)
        after-transform (assoc trimmed :intensities transformed-intensities)

        ;; Step 3: Smooth if requested
        smooth-fn (requiring-resolve 'maldi-clj.preprocessing/smooth-intensities)
        smoothed-intensities (smooth-fn (:intensities after-transform)
                                        smooth-method
                                        :half-window-size smooth-window)
        after-smooth (assoc after-transform :intensities smoothed-intensities)

        ;; Step 4: NEW - Baseline correction
        after-baseline (if (= baseline-method :none)
                         after-smooth
                         (remove-baseline after-smooth baseline-method
                                          :iterations baseline-iterations
                                          :half-window-size baseline-window))

        ;; Step 5: Calibrate (normalize)
        calibrate-fn (requiring-resolve 'maldi-clj.preprocessing/calibrate-intensities)
        calibrated-intensities (calibrate-fn (:intensities after-baseline) calibrate-method)
        final-spectrum (assoc after-baseline :intensities calibrated-intensities)]

    ;; Validate result
    (if (spectrum/valid-spectrum? final-spectrum)
      final-spectrum
      (throw (ex-info "Enhanced preprocessing with baseline correction produced invalid spectrum"
                      {:original spectrum :result final-spectrum})))))
