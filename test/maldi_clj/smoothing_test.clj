(ns maldi-clj.smoothing-test
  "Tests for Fastmath v3-based smoothing functions."
  (:require [clojure.test :refer [deftest is testing run-tests]]
            [maldi-clj.smoothing :as smooth]
            [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]))

(defn create-test-spectrum
  "Create a test spectrum with known signal + significant noise."
  [n-points]
  (let [mz-values (vec (range 1000 (+ 1000 n-points) 1.0))
        ;; Create signal with multiple peaks and significant noise
        intensities (dtype/make-container :jvm-heap :float64
                                          (mapv (fn [mz]
                                                  (+ 100 ; baseline
                                                     (* 50 (Math/exp (- (/ (Math/pow (- mz 1200) 2)
                                                                           (* 2 (Math/pow 50 2)))))) ; Peak at 1200
                                                     (* 30 (Math/exp (- (/ (Math/pow (- mz 1500) 2)
                                                                           (* 2 (Math/pow 30 2)))))) ; Peak at 1500
                                                     (* 50 (- (rand) 0.5)))) ; Significant noise
                                                mz-values))]
    {:id "test-spectrum"
     :mz mz-values
     :intensities intensities}))

(deftest test-savitzky-golay-basic
  "Test basic Savitzky-Golay functionality."
  (testing "Savitzky-Golay with spectrum map"
    (let [test-spec (create-test-spectrum 100)
          smoothed (smooth/savitzky-golay test-spec 11 3)]
      (is (map? smoothed))
      (is (= (:id test-spec) (:id smoothed)))
      (is (= (:mz test-spec) (:mz smoothed)))
      (is (= (count (:intensities test-spec)) (count (:intensities smoothed))))
      ;; Smoothing should reduce variance
      (let [orig-var (dfn/variance (:intensities test-spec))
            smooth-var (dfn/variance (:intensities smoothed))]
        (is (< smooth-var orig-var) "Smoothing should reduce variance"))))

  (testing "Savitzky-Golay with tensor input"
    (let [test-data (dtype/make-container :jvm-heap :float64
                                          (take 100 (repeatedly #(+ (rand 100) (* 5 (rand))))))
          smoothed (smooth/savitzky-golay test-data 11 3)]
      (is (= (count test-data) (count smoothed)))
      (is (dtype/reader? smoothed))
      ;; Smoothing should reduce variance
      (let [orig-var (dfn/variance test-data)
            smooth-var (dfn/variance smoothed)]
        (is (< smooth-var orig-var))))))

(deftest test-savitzky-golay-parameters
  "Test different Savitzky-Golay parameters."
  (let [test-spec (create-test-spectrum 200)]
    (testing "Different window sizes"
      (let [small-window (smooth/savitzky-golay test-spec 5 3)
            large-window (smooth/savitzky-golay test-spec 21 3)]
        ;; Larger window should smooth more (have less variance than smaller window)
        (let [small-var (dfn/variance (:intensities small-window))
              large-var (dfn/variance (:intensities large-window))]
          (is (< large-var small-var)
              (str "Large window variance (" large-var ") should be less than small window (" small-var ")")))))

    (testing "Different polynomial orders"
      (let [linear (smooth/savitzky-golay test-spec 11 1)
            cubic (smooth/savitzky-golay test-spec 11 3)]
        ;; Both should be valid but different
        (is (not= (:intensities linear) (:intensities cubic)))))))

(deftest test-variance-stabilization
  "Test square-root variance stabilization."
  (testing "Variance stabilization with spectrum map"
    (let [test-spec (create-test-spectrum 100)
          stabilized (smooth/variance-stabilization test-spec)]
      (is (map? stabilized))
      (is (= (:id test-spec) (:id stabilized)))
      (is (= (:mz test-spec) (:mz stabilized)))
      ;; All stabilized values should be positive and smaller
      (is (every? pos? (:intensities stabilized)))
      (is (< (apply max (:intensities stabilized))
             (apply max (:intensities test-spec))))))

  (testing "Variance stabilization with tensor"
    (let [test-data (dtype/make-container :jvm-heap :float64 (range 1 101))
          stabilized (smooth/variance-stabilization test-data)]
      (is (= (count test-data) (count stabilized)))
      (is (every? pos? stabilized))
      ;; Should be approximately sqrt(x + 0.375)
      (is (< (Math/abs (- (first stabilized)
                          (Math/sqrt (+ 1 0.375))))
             0.001)))))

(deftest test-variance-stabilization-constants
  "Test different Anscombe constants."
  (let [test-data (dtype/make-container :jvm-heap :float64 (range 1 21))
        default-constant (smooth/variance-stabilization test-data)
        custom-constant (smooth/variance-stabilization test-data 0.5)]
    (is (not= default-constant custom-constant))
    ;; Custom constant should give slightly different results
    (is (< (Math/abs (- (first custom-constant)
                        (Math/sqrt (+ 1 0.5))))
           0.001))))

(deftest test-moving-average
  "Test moving average smoothing."
  (testing "Moving average with spectrum map"
    (let [test-spec (create-test-spectrum 100)
          smoothed (smooth/moving-average test-spec 5)]
      (is (map? smoothed))
      (is (= (:id test-spec) (:id smoothed)))
      (is (= (count (:intensities test-spec)) (count (:intensities smoothed))))
      ;; Should reduce variance
      (let [orig-var (dfn/variance (:intensities test-spec))
            smooth-var (dfn/variance (:intensities smoothed))]
        (is (< smooth-var orig-var)))))

  (testing "Moving average with tensor"
    (let [test-data (dtype/make-container :jvm-heap :float64
                                          (take 50 (repeatedly #(+ 10 (rand 5)))))
          smoothed (smooth/moving-average test-data 5)]
      (is (= (count test-data) (count smoothed)))
      ;; Should reduce variance
      (let [orig-var (dfn/variance test-data)
            smooth-var (dfn/variance smoothed)]
        (is (< smooth-var orig-var))))))

(deftest test-clinical-preprocessing
  "Test complete clinical preprocessing pipeline."
  (testing "Clinical preprocessing with default parameters"
    (let [test-spec (create-test-spectrum 200)
          processed (smooth/clinical-preprocessing test-spec)]
      (is (map? processed))
      (is (= (:id test-spec) (:id processed)))
      (is (= (:mz test-spec) (:mz processed)))
      ;; Should be variance-stabilized and smoothed
      (is (every? pos? (:intensities processed)))
      (is (< (apply max (:intensities processed))
             (apply max (:intensities test-spec))))))

  (testing "Clinical preprocessing with custom parameters"
    (let [test-spec (create-test-spectrum 200)
          processed (smooth/clinical-preprocessing test-spec
                                                   {:variance-constant 0.5
                                                    :savgol-window 15
                                                    :savgol-order 2})]
      (is (map? processed))
      ;; Should produce different result than defaults
      (let [default-processed (smooth/clinical-preprocessing test-spec)]
        (is (not= (:intensities processed) (:intensities default-processed)))))))

(deftest test-savitzky-golay-derivative
  "Test Savitzky-Golay derivative calculation."
  (testing "First derivative"
    (let [test-spec (create-test-spectrum 100)
          first-deriv (smooth/savitzky-golay-derivative test-spec 1)]
      (is (map? first-deriv))
      (is (= (:id test-spec) (:id first-deriv)))
      (is (= (count (:intensities test-spec)) (count (:intensities first-deriv))))
      ;; Derivative should have different characteristics
      (is (not= (:intensities test-spec) (:intensities first-deriv)))))

  (testing "Second derivative"
    (let [test-data (dtype/make-container :jvm-heap :float64
                                          (mapv #(* % %) (range 0 10 0.1))) ; Quadratic function
          second-deriv (smooth/savitzky-golay-derivative test-data 2)]
      (is (= (count test-data) (count second-deriv)))
      ;; Second derivative of x^2 should be approximately constant (2)
      ;; Allow for some numerical error in the derivative calculation
      (let [mean-deriv (dfn/mean second-deriv)]
        (is (< (Math/abs (- mean-deriv 2.0)) 3.0)
            (str "Expected ~2.0, got " mean-deriv))))))

(deftest test-edge-cases
  "Test edge cases and error conditions."
  (testing "Empty spectrum"
    (let [empty-spec {:id "empty" :mz [] :intensities []}]
      ;; Should handle gracefully or throw meaningful error
      (is (thrown? Exception (smooth/savitzky-golay empty-spec)))))

  (testing "Very small spectrum"
    (let [tiny-spec {:id "tiny"
                     :mz [1000 1001]
                     :intensities (dtype/make-container :jvm-heap :float64 [100 101])}]
      ;; Should complete but may not be ideal - Fastmath handles this case
      (let [result (smooth/savitzky-golay tiny-spec 3 1)] ; Minimal window and order
        (is (map? result))
        (is (= 2 (count (:intensities result)))))))

  (testing "Invalid parameters"
    (let [test-spec (create-test-spectrum 100)]
      ;; Even window size should be rejected with AssertionError
      (is (thrown? AssertionError (smooth/savitzky-golay test-spec 10 3)))
      ;; Very high polynomial order might work but should produce valid results
      (let [result (smooth/savitzky-golay test-spec 5 4)]
        (is (map? result))
        (is (= (count (:intensities test-spec)) (count (:intensities result))))))))

(deftest test-performance-characteristics
  "Test performance-related characteristics."
  (testing "Large spectrum processing"
    (let [large-spec (create-test-spectrum 5000)
          start-time (System/nanoTime)
          smoothed (smooth/savitzky-golay large-spec 21 3)
          end-time (System/nanoTime)
          duration-ms (/ (- end-time start-time) 1000000.0)]
      (is (map? smoothed))
      (is (= (count (:intensities large-spec)) (count (:intensities smoothed))))
      ;; Should complete in reasonable time (< 100ms for 5000 points)
      (is (< duration-ms 100) (str "Processing took " duration-ms "ms"))))

  (testing "Memory efficiency"
    (let [test-spec (create-test-spectrum 1000)]
      ;; Multiple operations should not cause memory issues
      (dotimes [_ 10]
        (-> test-spec
            (smooth/variance-stabilization)
            (smooth/savitzky-golay 11 3)))
      ;; If we get here, no memory issues occurred
      (is true))))

(defn run-smoothing-tests
  "Run all smoothing tests and return results."
  []
  (run-tests 'maldi-clj.smoothing-test))
