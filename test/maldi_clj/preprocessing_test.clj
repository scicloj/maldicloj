(ns maldi-clj.preprocessing-test
  "Comprehensive tests for preprocessing pipeline functions."
  (:require [clojure.test :refer [deftest is testing run-tests]]
            [maldi-clj.preprocessing :as prep]
            [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.properties :as prop]
            [clojure.test.check.generators :as gen]))

;; =============================================================================
;; Test Data Generators
;; =============================================================================

(defn create-test-spectrum
  "Create a test spectrum with known characteristics."
  [n-points & {:keys [use-dtype? include-noise? baseline]
               :or {use-dtype? true include-noise? true baseline 10.0}}]
  (let [mz-values (vec (range 1000.0 (+ 1000.0 n-points) 1.0))
        intensities (mapv (fn [mz]
                            (+ baseline
                               (* 100 (Math/exp (- (/ (Math/pow (- mz 1500) 2) 50000))))
                               (if include-noise? (* 5 (- (rand) 0.5)) 0.0)))
                          mz-values)]
    (if use-dtype?
      (spectrum/create-spectrum "test-spectrum" mz-values intensities {})
      (spectrum/create-spectrum-legacy "test-spectrum" mz-values intensities {}))))

(defn create-noisy-intensities
  "Create noisy intensity data for testing."
  [n-points & {:keys [use-dtype?] :or {use-dtype? true}}]
  (let [intensities (mapv #(+ % (* 10 (rand))) (range 1 (inc n-points)))]
    (if use-dtype?
      (dtype/make-container :jvm-heap :float64 intensities)
      intensities)))

;; =============================================================================
;; Transformation Function Tests
;; =============================================================================

(deftest test-sqrt-stabilize
  "Test square root variance stabilization."
  (testing "Basic sqrt stabilization with dtype-next"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 4 9 16 25])
          result (prep/sqrt-stabilize test-data)]
      (is (dtype/reader? result))
      (is (every? identity
                  (map #(< (Math/abs (- %1 %2)) 0.001)
                       result [1.0 2.0 3.0 4.0 5.0])))))

  (testing "Sqrt stabilization with negative values"
    (let [test-data (dtype/make-container :jvm-heap :float64 [-1 0 1 4])
          result (prep/sqrt-stabilize test-data)]
      (is (every? #(>= % 0.0) result))
      (is (= 0.0 (first result))))) ; Negative becomes 0

  (testing "Sqrt stabilization with legacy vectors"
    (let [test-data [1 4 9 16 25]
          result (prep/sqrt-stabilize test-data)]
      (is (dtype/reader? result)) ; Returns dtype reader, not vector
      (is (every? identity
                  (map #(< (Math/abs (- %1 %2)) 0.001)
                       result [1.0 2.0 3.0 4.0 5.0]))))))

(deftest test-log-transforms
  "Test logarithmic transformations."
  (testing "Natural log transformation"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 (Math/E) (* (Math/E) (Math/E))])
          result (prep/log-transform test-data)]
      (is (< (Math/abs (- (nth result 0) 0.0)) 0.001))
      (is (< (Math/abs (- (nth result 1) 1.0)) 0.001))
      (is (< (Math/abs (- (nth result 2) 2.0)) 0.001))))

  (testing "Log10 transformation"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 10 100])
          result (prep/log10-transform test-data)]
      (is (< (Math/abs (- (nth result 0) 0.0)) 0.001))
      (is (< (Math/abs (- (nth result 1) 1.0)) 0.001))
      (is (< (Math/abs (- (nth result 2) 2.0)) 0.001))))

  (testing "Log2 transformation"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 2 4 8])
          result (prep/log2-transform test-data)]
      (is (< (Math/abs (- (nth result 0) 0.0)) 0.001))
      (is (< (Math/abs (- (nth result 1) 1.0)) 0.001))
      (is (< (Math/abs (- (nth result 2) 2.0)) 0.001))
      (is (< (Math/abs (- (nth result 3) 3.0)) 0.001))))

  (testing "Log transformation with zero/negative values"
    (let [test-data (dtype/make-container :jvm-heap :float64 [0 -1 1])
          result (prep/log-transform test-data)]
      (is (every? #(not (Double/isNaN %)) result))
      (is (every? #(not (Double/isInfinite %)) result)))))

(deftest test-transform-intensities
  "Test the general transform-intensities function."
  (let [test-data (dtype/make-container :jvm-heap :float64 [1 4 9 16])]
    (testing "Available transformation methods"
      (is (= test-data (prep/transform-intensities test-data :none)))
      (is (dtype/reader? (prep/transform-intensities test-data :sqrt)))
      (is (dtype/reader? (prep/transform-intensities test-data :log)))
      (is (dtype/reader? (prep/transform-intensities test-data :log10)))
      (is (dtype/reader? (prep/transform-intensities test-data :log2))))

    (testing "Invalid transformation method"
      (is (thrown? Exception (prep/transform-intensities test-data :invalid))))))

;; =============================================================================
;; Normalization Function Tests
;; =============================================================================

(deftest test-tic-normalize
  "Test Total Ion Current normalization."
  (testing "TIC normalization with dtype-next"
    (let [test-data (dtype/make-container :jvm-heap :float64 [10 20 30 40])
          result (prep/tic-normalize test-data)
          expected-tic 100.0]
      (is (dtype/reader? result))
      (is (< (Math/abs (- (dfn/sum result) 1.0)) 0.001))
      (is (every? identity
                  (map #(< (Math/abs (- %1 %2)) 0.001)
                       result [0.1 0.2 0.3 0.4])))))

  (testing "TIC normalization with zero sum"
    (let [test-data (dtype/make-container :jvm-heap :float64 [0 0 0])
          result (prep/tic-normalize test-data)]
      (is (= test-data result)))) ; Should return unchanged

  (testing "TIC normalization with legacy vectors"
    (let [test-data [10 20 30 40]
          result (prep/tic-normalize test-data)]
      (is (dtype/reader? result)) ; Returns dtype reader, not vector
      (is (< (Math/abs (- (reduce + result) 1.0)) 0.001)))))

(deftest test-median-normalize
  "Test median normalization."
  (testing "Median normalization with dtype-next"
    (let [test-data (dtype/make-container :jvm-heap :float64 [10 20 30 40])
          result (prep/median-normalize test-data)
          median-val 25.0] ; median of [10 20 30 40]
      (is (dtype/reader? result))
      (is (< (Math/abs (- (dfn/median result) 1.0)) 0.1)))) ; Allow some tolerance

  (testing "Median normalization with zero median"
    (let [test-data (dtype/make-container :jvm-heap :float64 [0 0 0])
          result (prep/median-normalize test-data)]
      (is (= test-data result))))

  (testing "Empty data handling"
    (let [test-data (dtype/make-container :jvm-heap :float64 [])
          result (prep/median-normalize test-data)]
      (is (= 0 (dtype/ecount result))))))

(deftest test-calibrate-intensities
  "Test the general calibrate-intensities function."
  (let [test-data (dtype/make-container :jvm-heap :float64 [10 20 30 40])]
    (testing "Available calibration methods"
      (is (= test-data (prep/calibrate-intensities test-data :none)))
      (is (dtype/reader? (prep/calibrate-intensities test-data :tic)))
      (is (dtype/reader? (prep/calibrate-intensities test-data :median))))

    (testing "Invalid calibration method"
      (is (thrown? Exception (prep/calibrate-intensities test-data :invalid))))))

;; =============================================================================
;; Smoothing Function Tests
;; =============================================================================

(deftest test-moving-average-smooth
  "Test moving average smoothing."
  (testing "Simple moving average"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 2 3 4 5 4 3 2 1])
          result (prep/moving-average-smooth test-data 1)]
      (is (dtype/reader? result))
      (is (= (dtype/ecount test-data) (dtype/ecount result)))
;; Middle value should be average of [4,5,4] = 4.33...
      (is (< (Math/abs (- (nth result 4) 4.333333333333333)) 0.001))))

  (testing "Weighted moving average"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 2 3 4 5])
          result (prep/moving-average-smooth test-data 2 :weighted? true)]
      (is (dtype/reader? result))
      (is (= (dtype/ecount test-data) (dtype/ecount result)))))

  (testing "Edge case: small window"
    (let [test-data (dtype/make-container :jvm-heap :float64 [1 2])
          result (prep/moving-average-smooth test-data 5)] ; Window larger than data
      (is (= (dtype/ecount test-data) (dtype/ecount result)))))

  (testing "Legacy vector support"
    (let [test-data [1 2 3 4 5]
          result (prep/moving-average-smooth test-data 1)]
      (is (dtype/reader? result)) ; Returns dtype reader, not vector
      (is (= (count test-data) (dtype/ecount result))))))

(deftest test-smooth-intensities
  "Test the general smooth-intensities function."
  (let [test-data (dtype/make-container :jvm-heap :float64 [1 5 2 6 3 7 4])]
    (testing "Available smoothing methods"
      (is (= test-data (prep/smooth-intensities test-data :none)))
      (is (dtype/reader? (prep/smooth-intensities test-data :moving-average)))
      (is (dtype/reader? (prep/smooth-intensities test-data :weighted-moving-average))))

    (testing "Custom window size"
      (let [result1 (prep/smooth-intensities test-data :moving-average :half-window-size 1)
            result2 (prep/smooth-intensities test-data :moving-average :half-window-size 3)]
        (is (not= result1 result2))))

    (testing "Invalid smoothing method"
      (is (thrown? Exception (prep/smooth-intensities test-data :invalid))))))

;; =============================================================================
;; Spectrum Manipulation Tests
;; =============================================================================

(deftest test-trim-spectrum
  "Test spectrum trimming functionality."
  (testing "Basic trimming with dtype spectrum"
    (let [test-spectrum (create-test-spectrum 1000 :use-dtype? true)
          trimmed (prep/trim-spectrum test-spectrum 1200.0 1800.0)
          trimmed-mz (:mz-values trimmed)]
      (is (spectrum/dtype-spectrum? trimmed))
      (is (every? #(<= 1200.0 % 1800.0) trimmed-mz))
      (is (< (dtype/ecount trimmed-mz) (dtype/ecount (:mz-values test-spectrum))))))

  (testing "Trimming with legacy spectrum"
    (let [test-spectrum (create-test-spectrum 100 :use-dtype? false)
          trimmed (prep/trim-spectrum test-spectrum 1020.0 1080.0)]
      (is (map? trimmed)) ; Just check it's a valid map
      (is (contains? trimmed :mz-values))
      (is (contains? trimmed :intensities))
      (is (every? #(<= 1020.0 % 1080.0) (:mz-values trimmed)))))

  (testing "Trimming to empty range"
    (let [test-spectrum (create-test-spectrum 100)
          trimmed (prep/trim-spectrum test-spectrum 5000.0 6000.0)] ; Outside range
      (is (zero? (dtype/ecount (:mz-values trimmed))))
      (is (zero? (dtype/ecount (:intensities trimmed))))))

  (testing "No trimming needed"
    (let [test-spectrum (create-test-spectrum 100)
          trimmed (prep/trim-spectrum test-spectrum 0.0 10000.0)]
      (is (= (dtype/ecount (:mz-values test-spectrum))
             (dtype/ecount (:mz-values trimmed)))))))

(deftest test-total-ion-current
  "Test TIC calculation."
  (testing "TIC with dtype spectrum"
    (let [test-spectrum (create-test-spectrum 100 :use-dtype? true)
          tic (prep/total-ion-current test-spectrum)]
      (is (number? tic))
      (is (pos? tic))))

  (testing "TIC with legacy spectrum"
    (let [test-spectrum (create-test-spectrum 100 :use-dtype? false)
          tic (prep/total-ion-current test-spectrum)]
      (is (number? tic))
      (is (pos? tic))))

  (testing "TIC with zero intensities"
    (let [test-spectrum (spectrum/create-spectrum "zero" [1000 1001] [0 0] {})]
      (is (zero? (prep/total-ion-current test-spectrum))))))

(deftest test-is-empty-spectrum
  "Test empty spectrum detection."
  (testing "Non-empty spectrum"
    (let [test-spectrum (create-test-spectrum 100)]
      (is (not (prep/is-empty-spectrum? test-spectrum)))))

  (testing "Empty spectrum"
    (let [empty-spectrum (spectrum/create-spectrum "empty" [] [] {})]
      (is (prep/is-empty-spectrum? empty-spectrum)))))

;; =============================================================================
;; Integration Pipeline Tests
;; =============================================================================

(deftest test-preprocess-spectrum-enhanced
  "Test enhanced preprocessing pipeline."
  (testing "Complete preprocessing pipeline"
    (let [test-spectrum (create-test-spectrum 200)
          result (prep/preprocess-spectrum-enhanced test-spectrum
                                                    :transform-method :sqrt
                                                    :smooth-method :moving-average
                                                    :calibrate-method :tic
                                                    :smooth-window 5
                                                    :trim-range [1100.0 1900.0])]
      (is (spectrum/valid-spectrum? result))
      (is (< (dtype/ecount (:mz-values result))
             (dtype/ecount (:mz-values test-spectrum)))) ; Should be trimmed
      (is (< (Math/abs (- (dfn/sum (:intensities result)) 1.0)) 0.001)))) ; Should be TIC normalized

  (testing "No operations pipeline"
    (let [test-spectrum (create-test-spectrum 100)
          result (prep/preprocess-spectrum-enhanced test-spectrum)]
      (is (spectrum/valid-spectrum? result))
      ;; Should be essentially unchanged except for potential dtype conversion
      (is (= (dtype/ecount (:mz-values test-spectrum))
             (dtype/ecount (:mz-values result))))))

  (testing "Invalid parameters"
    (let [test-spectrum (create-test-spectrum 100)]
      (is (thrown? Exception
                   (prep/preprocess-spectrum-enhanced test-spectrum
                                                      :transform-method :invalid)))
      (is (thrown? Exception
                   (prep/preprocess-spectrum-enhanced test-spectrum
                                                      :smooth-method :invalid)))))

  (testing "Legacy spectrum conversion"
    (let [legacy-spectrum (create-test-spectrum 100 :use-dtype? false)
          result (prep/preprocess-spectrum-enhanced legacy-spectrum
                                                    :force-dtype-conversion? true)]
      (is (spectrum/dtype-spectrum? result)))))

(deftest test-preprocess-spectrum-legacy
  "Test legacy preprocessing function."
  (testing "Legacy API compatibility"
    (let [test-spectrum (create-test-spectrum 100)
          result (prep/preprocess-spectrum test-spectrum
                                           :sqrt-stabilize? true
                                           :tic-normalize? true)]
      (is (spectrum/valid-spectrum? result))
      ;; Check that sqrt and TIC were applied
      (is (< (Math/abs (- (dfn/sum (:intensities result)) 1.0)) 0.001))))

  (testing "No operations"
    (let [test-spectrum (create-test-spectrum 100)
          result (prep/preprocess-spectrum test-spectrum)]
      (is (spectrum/valid-spectrum? result))))

  (testing "With legacy input spectrum"
    (let [legacy-spectrum (create-test-spectrum 100 :use-dtype? false)
          result (prep/preprocess-spectrum legacy-spectrum
                                           :sqrt-stabilize? true)]
      (is (spectrum/valid-spectrum? result)))))

;; =============================================================================
;; Performance and Edge Case Tests
;; =============================================================================

(deftest test-performance-characteristics
  "Test performance-related characteristics."
  (testing "Large spectrum processing"
    (let [large-spectrum (create-test-spectrum 5000)
          start-time (System/nanoTime)
          result (prep/preprocess-spectrum-enhanced large-spectrum
                                                    :transform-method :sqrt
                                                    :smooth-method :moving-average
                                                    :calibrate-method :tic)
          end-time (System/nanoTime)
          duration-ms (/ (- end-time start-time) 1000000.0)]
      (is (spectrum/valid-spectrum? result))
      ;; Should complete in reasonable time
      (is (< duration-ms 1000) (str "Processing took " duration-ms "ms"))))

  (testing "Memory efficiency"
    (let [test-spectrum (create-test-spectrum 1000)]
      ;; Multiple operations should not cause memory issues
      (dotimes [_ 5]
        (prep/preprocess-spectrum-enhanced test-spectrum
                                           :transform-method :sqrt
                                           :smooth-method :moving-average
                                           :calibrate-method :tic))
      ;; If we get here, no memory issues occurred
      (is true))))

(deftest test-edge-cases
  "Test edge cases and error conditions."
  (testing "Empty spectrum handling"
    (let [empty-spectrum (spectrum/create-spectrum "empty" [] [] {})]
      (is (spectrum/valid-spectrum?
           (prep/preprocess-spectrum-enhanced empty-spectrum)))))

  (testing "Single point spectrum"
    (let [tiny-spectrum (spectrum/create-spectrum "tiny" [1000.0] [100.0] {})]
      (is (spectrum/valid-spectrum?
           (prep/preprocess-spectrum-enhanced tiny-spectrum
                                              :smooth-method :moving-average)))))

  (testing "All zero intensities"
    (let [zero-spectrum (spectrum/create-spectrum "zeros" [1000 1001 1002] [0 0 0] {})]
      (is (spectrum/valid-spectrum?
           (prep/preprocess-spectrum-enhanced zero-spectrum
                                              :calibrate-method :tic)))))

  (testing "Very large values"
    (let [large-spectrum (spectrum/create-spectrum "large" [1000 1001] [1e10 1e10] {})]
      (is (spectrum/valid-spectrum?
           (prep/preprocess-spectrum-enhanced large-spectrum
                                              :transform-method :log))))))

;; =============================================================================
;; Property-Based Tests
;; =============================================================================

(defspec test-transform-preserves-length 50
  (prop/for-all [intensities (gen/vector (gen/double* {:min 0.1 :max 1000.0}) 10 100)]
                (let [input (dtype/make-container :jvm-heap :float64 intensities)
                      sqrt-result (prep/sqrt-stabilize input)
                      log-result (prep/log-transform input)]
                  (and (= (dtype/ecount input) (dtype/ecount sqrt-result))
                       (= (dtype/ecount input) (dtype/ecount log-result))))))

(defspec test-normalization-properties 50
  (prop/for-all [intensities (gen/vector (gen/double* {:min 1.0 :max 1000.0}) 5 50)]
                (let [input (dtype/make-container :jvm-heap :float64 intensities)
                      tic-result (prep/tic-normalize input)]
                  (and (= (dtype/ecount input) (dtype/ecount tic-result))
                       (< (Math/abs (- (dfn/sum tic-result) 1.0)) 0.001)))))

(defspec test-preprocessing-pipeline-validity 30
  (prop/for-all [n-points (gen/choose 10 100)
                 transform-method (gen/elements [:none :sqrt :log])
                 calibrate-method (gen/elements [:none :tic])]
                (let [test-spectrum (create-test-spectrum n-points)
                      result (prep/preprocess-spectrum-enhanced test-spectrum
                                                                :transform-method transform-method
                                                                :calibrate-method calibrate-method)]
                  (spectrum/valid-spectrum? result))))

;; =============================================================================
;; Test Runner
;; =============================================================================

(defn run-preprocessing-tests
  "Run all preprocessing tests and return results."
  []
  (run-tests 'maldi-clj.preprocessing-test))
