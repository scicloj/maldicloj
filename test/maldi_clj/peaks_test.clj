(ns maldi-clj.peaks-test
  (:require [clojure.test :refer :all]
            [maldi-clj.peaks :as peaks]))

(deftest test-basic-peak-detection
  (testing "Basic peak detection functionality"
    (let [mz-values [100.0 101.0 102.0 103.0 104.0]
          intensities [10.0 50.0 30.0 80.0 15.0]
          peaks (peaks/detect-peaks mz-values intensities {:snr-threshold 2.0})]
      
      (is (not (empty? peaks)))
      (is (every? #(contains? % :mz) peaks))
      (is (every? #(contains? % :intensity) peaks))
      (is (every? #(contains? % :index) peaks))
      (is (every? #(contains? % :snr) peaks)))))

(deftest test-optimization-consistency
  (testing "Optimized vs standard implementation consistency"
    (let [mz-values (vec (range 0 2000 0.1))
          intensities (vec (repeatedly 20000 #(+ 10 (rand 1000))))
          
          standard-peaks (peaks/detect-peaks mz-values intensities 
                                            {:use-optimizations? false
                                             :snr-threshold 3.0})
          optimized-peaks (peaks/detect-peaks mz-values intensities 
                                             {:use-optimizations? true
                                              :snr-threshold 3.0})]
      
      (is (= (count standard-peaks) (count optimized-peaks)))
      (is (= (set (map :index standard-peaks)) 
             (set (map :index optimized-peaks)))))))

(deftest test-performance-benchmark
  (testing "Performance benchmarking functionality"
    (let [mz-values (vec (range 1000))
          intensities (vec (repeatedly 1000 #(+ 10 (rand 1000))))
          benchmark (peaks/benchmark-peak-detection mz-values intensities {})]
      
      (is (contains? benchmark :standard))
      (is (contains? benchmark :optimized))
      (is (contains? benchmark :speedup))
      (is (pos? (:speedup benchmark)))
      (is (pos? (get-in benchmark [:standard :time-ms])))
      (is (pos? (get-in benchmark [:optimized :time-ms]))))))

(deftest test-batch-processing
  (testing "Batch peak detection"
    (let [spectra-map {"spec1" {:mz [100 101 102 103 104]
                               :intensities [10 50 30 80 15]}
                       "spec2" {:mz [200 201 202 203 204] 
                               :intensities [20 60 40 90 25]}}
          results (peaks/detect-peaks-batch spectra-map {})]
      
      (is (= (set (keys results)) #{"spec1" "spec2"}))
      (is (every? vector? (vals results)))
      (is (every? #(every? map? %) (vals results))))))

(deftest test-parallel-batch-processing
  (testing "Parallel batch processing"
    (let [large-spectra-map (into {} 
                                  (for [i (range 150)]
                                    [(str "spec" i) 
                                     {:mz (vec (range 100))
                                      :intensities (vec (repeatedly 100 #(rand 1000)))}]))
          
          sequential-results (peaks/detect-peaks-batch large-spectra-map 
                                                      {:parallel? false})
          parallel-results (peaks/detect-peaks-batch large-spectra-map 
                                                    {:parallel? true
                                                     :chunk-size 50})]
      
      (is (= (count sequential-results) (count parallel-results)))
      (is (= (set (keys sequential-results)) (set (keys parallel-results)))))))

(deftest test-centroiding
  (testing "Peak centroiding functionality"
    (let [mz-values [99.9 100.0 100.1 100.2]
          intensities [30.0 100.0 80.0 20.0]
          peaks [{:mz 100.0 :intensity 100.0 :index 1}]
          centroided (peaks/centroid-peaks mz-values intensities peaks {})]
      
      (is (= (count centroided) 1))
      (is (contains? (first centroided) :centroided))
      (is (:centroided (first centroided)))
      ;; Centroided m/z should be weighted average, so between 100.0 and 100.1
      (is (and (>= (:mz (first centroided)) 100.0)
               (<= (:mz (first centroided)) 100.1))))))

(deftest test-memory-estimation
  (testing "Memory usage estimation"
    (let [estimation (peaks/estimate-memory-usage 10000)]
      
      (is (contains? estimation :base-mb))
      (is (contains? estimation :temp-mb)) 
      (is (contains? estimation :result-mb))
      (is (contains? estimation :total-mb))
      (is (pos? (:total-mb estimation))))))

(deftest test-performance-profile  
  (testing "Performance profiling across data sizes"
    (let [profile (peaks/performance-profile :sizes [100 500])]
      
      (is (= (count profile) 2))
      (is (every? #(contains? % :size) profile))
      (is (every? #(contains? % :speedup) profile))
      (is (= (map :size profile) [100 500])))))

(deftest test-edge-cases
  (testing "Edge cases and error handling"
    
    ;; Empty data
    (is (empty? (peaks/detect-peaks [] [])))
    
    ;; Mismatched lengths
    (is (thrown? Exception 
                (peaks/detect-peaks [1 2 3] [1 2])))
    
    ;; Single data point
    (let [peaks (peaks/detect-peaks [100.0] [50.0])]
      (is (<= (count peaks) 1)))
    
    ;; All zeros
    (let [peaks (peaks/detect-peaks [100 101 102] [0 0 0])]
      (is (empty? peaks)))))

(deftest test-filtering
  (testing "Peak filtering functionality"
    (let [peaks [{:mz 100.0 :intensity 50.0 :snr 2.5}
                 {:mz 101.0 :intensity 100.0 :snr 5.0}
                 {:mz 102.0 :intensity 25.0 :snr 1.5}]
          
          filtered-snr (peaks/filter-peaks peaks {:min-snr 3.0})
          filtered-intensity (peaks/filter-peaks peaks {:intensity-threshold 40.0})
          filtered-count (peaks/filter-peaks peaks {:max-peaks 2})]
      
      (is (= (count filtered-snr) 1))
      (is (= (:mz (first filtered-snr)) 101.0))
      
      (is (= (count filtered-intensity) 2))
      
      (is (<= (count filtered-count) 2)))))

(run-tests)(ns maldi-clj.peaks-test
  "Tests for peak detection algorithms"
  (:require [clojure.test :refer [deftest is testing]]
            [maldi-clj.peaks :as peaks]
            [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]))

;; ==============================================================================
;; TEST UTILITIES
;; ==============================================================================

(defn create-synthetic-spectrum-with-peaks
  "Create synthetic spectrum with known peaks for testing"
  [n-points]
  (let [mz-values (mapv #(+ 100 (* % 0.1)) (range n-points))
        intensities (mapv (fn [i]
                            (+ 50 ; base intensity
                               (* 10 (rand)) ; random noise
                               (cond
                                 ;; Peak at index 200 (m/z ~120)
                                 (< (Math/abs (- i 200)) 3) (* 1000 (Math/exp (- (/ (* (- i 200) (- i 200)) 4))))
                                 ;; Peak at index 500 (m/z ~150) 
                                 (< (Math/abs (- i 500)) 2) (* 800 (Math/exp (- (/ (* (- i 500) (- i 500)) 2))))
                                 ;; Peak at index 700 (m/z ~170)
                                 (< (Math/abs (- i 700)) 4) (* 600 (Math/exp (- (/ (* (- i 700) (- i 700)) 8))))
                                 :else 0)))
                          (range n-points))]
    (spectrum/create-spectrum "synthetic-peaks" mz-values intensities
                              {:synthetic true
                               :expected-peaks [{:mz 120.0 :index 200}
                                                {:mz 150.0 :index 500}
                                                {:mz 170.0 :index 700}]})))

(defn create-noise-spectrum
  "Create spectrum with only noise for noise estimation testing"
  [n-points noise-level]
  (let [mz-values (mapv #(+ 100 (* % 0.1)) (range n-points))
        intensities (mapv (fn [_] (+ noise-level (* noise-level 0.2 (- (rand) 0.5))))
                          (range n-points))]
    (spectrum/create-spectrum "noise-only" mz-values intensities {:noise-level noise-level})))

;; ==============================================================================
;; NOISE ESTIMATION TESTS
;; ==============================================================================

(deftest test-noise-estimation-mad
  (testing "MAD noise estimation"
    (let [;; Test with known noise level
          noise-spectrum (create-noise-spectrum 1000 100)
          estimated-noise (peaks/estimate-noise-mad (:intensities noise-spectrum))]

      (testing "MAD should estimate noise reasonably"
        (is (> estimated-noise 0) "Noise estimate should be positive")
        (is (< 5 estimated-noise 50) "Noise estimate should be reasonable for synthetic data"))

      (testing "MAD should work with dtype containers"
        (is (number? estimated-noise) "Should return a single number"))

      (testing "MAD should work with legacy vectors"
        (let [legacy-intensities (dtype/->vector (:intensities noise-spectrum))
              legacy-noise (peaks/estimate-noise-mad legacy-intensities)]
          (is (< (Math/abs (- estimated-noise legacy-noise)) 1.0)
              "dtype and legacy implementations should give similar results"))))))

(deftest test-noise-estimation-local-mad
  (testing "Local MAD noise estimation"
    (let [noise-spectrum (create-noise-spectrum 500 100)
          local-noise (peaks/estimate-noise-local-mad (:intensities noise-spectrum)
                                                      :window-half-size 25)]

      (testing "Local MAD should return vector of estimates"
        (is (dtype/reader? local-noise) "Should return dtype container")
        (is (= (spectrum/spectrum-length noise-spectrum) (dtype/ecount local-noise))
            "Should have noise estimate for each point"))

      (testing "Local noise estimates should be reasonable"
        (let [noise-values (dtype/->vector local-noise)]
          (is (every? #(> % 0) noise-values) "All noise estimates should be positive")
          (is (< 5 (dfn/mean local-noise) 50) "Mean noise should be reasonable"))))))

;; ==============================================================================
;; LOCAL MAXIMA DETECTION TESTS
;; ==============================================================================

(deftest test-local-maxima-detection
  (testing "Local maxima detection"
    (let [;; Simple test data with clear maxima
          intensities (dtype/make-container :jvm-heap :float32
                                            [1, 2, 5, 8, 5, 2, 1, 3, 9, 12, 9, 3, 1, 2, 1])
          maxima (peaks/find-local-maxima intensities :window-half-size 1)]

      (testing "Should find correct local maxima"
        (is (= [3 9 13] maxima) "Should find peaks at indices 3, 9, and 13"))

      (testing "Should respect minimum intensity threshold"
        (let [filtered-maxima (peaks/find-local-maxima intensities
                                                       :window-half-size 1
                                                       :min-intensity 10)]
          (is (= [9] filtered-maxima) "Should only find peak above intensity 10")))

      (testing "Should work with different window sizes"
        (let [large-window-maxima (peaks/find-local-maxima intensities :window-half-size 3)]
          (is (< (count large-window-maxima) (count maxima)) "Larger window should find fewer maxima")
          (is (some #(= % 9) large-window-maxima) "Should still find the highest peak"))))))

(deftest test-local-maxima-with-synthetic-peaks
  (testing "Local maxima on synthetic peak data"
    (let [peak-spectrum (create-synthetic-spectrum-with-peaks 1000)
          intensities (:intensities peak-spectrum)
          maxima (peaks/find-local-maxima intensities
                                          :window-half-size 5
                                          :min-intensity 100)]

      (testing "Should find peaks near expected locations"
        ;; Expected peaks around indices 200, 500, 700
        (let [peak-regions {200 [190 210] 500 [490 510] 700 [690 710]}
              found-in-regions (reduce-kv
                                (fn [acc expected-idx [min-idx max-idx]]
                                  (assoc acc expected-idx
                                         (some #(<= min-idx % max-idx) maxima)))
                                {}
                                peak-regions)]
          (is (every? identity (vals found-in-regions))
              "Should find maxima in all expected peak regions"))))))

;; ==============================================================================
;; SNR CALCULATION TESTS
;; ==============================================================================

(deftest test-snr-calculation
  (testing "SNR calculation"
    (testing "Basic SNR calculation"
      (is (= 3.0 (peaks/calculate-snr 100 25)) "SNR = (100-25)/25 = 3")
      (is (= 9.0 (peaks/calculate-snr 100 10)) "SNR = (100-10)/10 = 9"))

    (testing "Edge cases"
      (is (= Double/POSITIVE_INFINITY (peaks/calculate-snr 100 0))
          "Zero noise should give infinite SNR")
      (is (= -0.5 (peaks/calculate-snr 50 100))
          "Signal below noise should give negative SNR"))

    (testing "SNR vector calculation"
      (let [intensities (dtype/make-container :jvm-heap :float32 [100 200 50])
            noise 50
            snr-values (peaks/calculate-snr-vector intensities noise)]
        (is (= [1.0 3.0 0.0] (dtype/->vector snr-values))
            "Should calculate SNR for each intensity value")))))

;; ==============================================================================
;; FULL PEAK DETECTION TESTS
;; ==============================================================================

(deftest test-peak-detection-basic
  (testing "Basic peak detection functionality"
    (let [peak-spectrum (create-synthetic-spectrum-with-peaks 1000)
          detected-peaks (peaks/detect-peaks peak-spectrum
                                             :snr-threshold 10.0
                                             :window-half-size 5
                                             :baseline-corrected? true)]

      (testing "Should return valid MassPeaks structure"
        (is (map? detected-peaks) "Should return a map")
        (is (string? (:id detected-peaks)) "Should have string ID")
        (is (vector? (:peaks detected-peaks)) "Should have vector of peaks")
        (is (map? (:metadata detected-peaks)) "Should have metadata map"))

      (testing "Should detect expected peaks"
        (is (> (peaks/peak-count detected-peaks) 0) "Should detect some peaks")
        (let [major-peaks (filter #(> (:snr %) 50) (:peaks detected-peaks))]
          (is (>= (count major-peaks) 2) "Should detect at least 2 major peaks")))

      (testing "Peak structure should be valid"
        (let [first-peak (first (:peaks detected-peaks))]
          (is (number? (:mz first-peak)) "Peak should have m/z value")
          (is (number? (:intensity first-peak)) "Peak should have intensity")
          (is (number? (:snr first-peak)) "Peak should have SNR")
          (is (number? (:noise first-peak)) "Peak should have noise estimate")
          (is (int? (:index first-peak)) "Peak should have index"))))))

(deftest test-peak-detection-parameters
  (testing "Peak detection parameter effects"
    (let [peak-spectrum (create-synthetic-spectrum-with-peaks 1000)]

      (testing "SNR threshold effects"
        (let [low-snr-peaks (peaks/detect-peaks peak-spectrum :snr-threshold 1.0
                                                :baseline-corrected? true)
              high-snr-peaks (peaks/detect-peaks peak-spectrum :snr-threshold 20.0
                                                 :baseline-corrected? true)]
          (is (> (peaks/peak-count low-snr-peaks) (peaks/peak-count high-snr-peaks))
              "Lower SNR threshold should detect more peaks")))

      (testing "Window size effects"
        (let [small-window-peaks (peaks/detect-peaks peak-spectrum :window-half-size 2
                                                     :baseline-corrected? true)
              large-window-peaks (peaks/detect-peaks peak-spectrum :window-half-size 20
                                                     :baseline-corrected? true)]
          (is (>= (peaks/peak-count small-window-peaks) (peaks/peak-count large-window-peaks))
              "Smaller window should detect same or more peaks")))

      (testing "Noise method effects"
        (let [mad-peaks (peaks/detect-peaks peak-spectrum :noise-method :mad
                                            :baseline-corrected? true)
              local-mad-peaks (peaks/detect-peaks peak-spectrum :noise-method :local-mad
                                                  :baseline-corrected? true)]
          (is (> (peaks/peak-count mad-peaks) 0) "MAD method should detect peaks")
          (is (> (peaks/peak-count local-mad-peaks) 0) "Local MAD method should detect peaks"))))))

;; ==============================================================================
;; PEAK UTILITIES TESTS
;; ==============================================================================

(deftest test-peak-utilities
  (testing "Peak utility functions"
    (let [peak-spectrum (create-synthetic-spectrum-with-peaks 500)
          detected-peaks (peaks/detect-peaks peak-spectrum :snr-threshold 5.0
                                             :baseline-corrected? true)]

      (testing "Peak count"
        (is (= (count (:peaks detected-peaks)) (peaks/peak-count detected-peaks))
            "Peak count should match vector length"))

      (testing "Peak value extraction"
        (let [intensities (peaks/get-peak-intensities detected-peaks)
              mz-values (peaks/get-peak-mz-values detected-peaks)
              snr-values (peaks/get-peak-snr-values detected-peaks)]
          (is (= (count intensities) (peaks/peak-count detected-peaks))
              "Should extract correct number of intensities")
          (is (every? number? intensities) "All intensities should be numbers")
          (is (every? number? mz-values) "All m/z values should be numbers")
          (is (every? number? snr-values) "All SNR values should be numbers")))

      (testing "Peak filtering"
        (let [high-snr-peaks (peaks/filter-peaks-by-snr detected-peaks 20.0)
              high-intensity-peaks (peaks/filter-peaks-by-intensity detected-peaks 500)]
          (is (<= (peaks/peak-count high-snr-peaks) (peaks/peak-count detected-peaks))
              "Filtering should reduce or maintain peak count")
          (is (every? #(>= (:snr %) 20.0) (:peaks high-snr-peaks))
              "Filtered peaks should meet SNR criterion")
          (is (every? #(>= (:intensity %) 500) (:peaks high-intensity-peaks))
              "Filtered peaks should meet intensity criterion")))

      (testing "Peak statistics"
        (let [stats (peaks/peak-statistics detected-peaks)]
          (is (map? stats) "Should return statistics map")
          (is (number? (:count stats)) "Should have peak count")
          (is (map? (:intensity-stats stats)) "Should have intensity statistics")
          (is (map? (:snr-stats stats)) "Should have SNR statistics"))))))

;; ==============================================================================
;; EDGE CASES AND ERROR HANDLING
;; ==============================================================================

(deftest test-edge-cases
  (testing "Edge cases and error handling"
    (testing "Empty spectrum"
      (let [empty-spectrum (spectrum/create-spectrum "empty" [] [] {})]
        (is (thrown? Exception (peaks/detect-peaks empty-spectrum))
            "Should handle empty spectrum gracefully")))

    (testing "Single point spectrum"
      (let [single-point (spectrum/create-spectrum "single" [100.0] [1000.0] {})]
        (let [detected-peaks (peaks/detect-peaks single-point :baseline-corrected? true)]
          (is (>= (peaks/peak-count detected-peaks) 0) "Should handle single point"))))

    (testing "Invalid noise method"
      (let [test-spectrum (create-synthetic-spectrum-with-peaks 100)]
        (is (thrown? Exception
                     (peaks/detect-peaks test-spectrum :noise-method :invalid))
            "Should throw exception for invalid noise method")))

    (testing "Noise estimation with constant values"
      (let [constant-intensities (dtype/make-container :jvm-heap :float32
                                                       (repeat 100 50.0))
            noise (peaks/estimate-noise-mad constant-intensities)]
        (is (>= noise 0) "Should handle constant values without error")))))

;; ==============================================================================
;; PERFORMANCE TESTS
;; ==============================================================================

(deftest test-performance-basic
  (testing "Basic performance characteristics"
    (let [large-spectrum (create-synthetic-spectrum-with-peaks 10000)]

      (testing "Large spectrum processing"
        (let [start-time (System/nanoTime)
              detected-peaks (peaks/detect-peaks large-spectrum :baseline-corrected? true)
              end-time (System/nanoTime)
              duration-ms (/ (- end-time start-time) 1000000.0)]

          (is (> (peaks/peak-count detected-peaks) 0) "Should detect peaks in large spectrum")
          (is (< duration-ms 5000) "Should process 10k points in reasonable time")
          (println (format "Peak detection on 10k points took %.2f ms" duration-ms)))))))

;; ==============================================================================
;; CONVENIENCE FUNCTION TESTS
;; ==============================================================================

(deftest test-convenience-functions
  (testing "Convenience wrapper functions"
    (let [test-spectrum (create-synthetic-spectrum-with-peaks 500)]

      (testing "estimate-noise convenience function"
        (let [mad-noise (peaks/estimate-noise test-spectrum :mad)
              local-mad-noise (peaks/estimate-noise test-spectrum :local-mad :window-half-size 25)]
          (is (number? mad-noise) "MAD noise should return single number")
          (is (dtype/reader? local-mad-noise) "Local MAD should return container")
          (is (= (spectrum/spectrum-length test-spectrum) (dtype/ecount local-mad-noise))
              "Local noise should have estimate for each point")))

      (testing "Invalid noise method"
        (is (thrown? Exception (peaks/estimate-noise test-spectrum :invalid))
            "Should throw exception for invalid method")))))

;; ==============================================================================
;; TEST RUNNER
;; ==============================================================================

(defn run-peak-tests
  "Run all peak detection tests"
  []
  (clojure.test/run-tests 'maldi-clj.peaks-test))

(defn run-peak-performance-tests
  "Run performance-focused tests"
  []
  (test-performance-basic))
