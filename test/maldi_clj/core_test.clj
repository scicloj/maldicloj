(ns maldi-clj.core-test
  "Comprehensive tests for the core API functions."
  (:require [clojure.test :refer [deftest is testing run-tests]]
            [maldi-clj.core :as core]
            [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [tech.v3.dataset :as ds]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.properties :as prop]
            [clojure.test.check.generators :as gen]))

;; =============================================================================
;; Test Data Generators
;; =============================================================================

(defn create-test-data
  "Create test m/z and intensity data."
  [n-points & {:keys [with-peaks? baseline noise-level]
               :or {with-peaks? true baseline 50.0 noise-level 5.0}}]
  (let [mz-values (vec (range 1000.0 (+ 1000.0 n-points) 1.0))
        intensities (mapv (fn [mz]
                            (+ baseline
                               (if with-peaks?
                                 (+ (* 100 (Math/exp (- (/ (Math/pow (- mz 1200) 2) 10000))))
                                    (* 80 (Math/exp (- (/ (Math/pow (- mz 1500) 2) 8000))))
                                    (* 60 (Math/exp (- (/ (Math/pow (- mz 1800) 2) 12000)))))
                                 0.0)
                               (* noise-level (- (rand) 0.5))))
                          mz-values)]
    [mz-values intensities]))

;; =============================================================================
;; Core Spectrum Creation Tests
;; =============================================================================

(deftest test-create-spectrum
  "Test core spectrum creation function."
  (testing "Basic spectrum creation with defaults"
    (let [[mz-vals intensities] (create-test-data 100)
          spectrum (core/create-spectrum "test-001" mz-vals intensities
                                         {:sample-type "test"})]
      (is (spectrum/valid-spectrum? spectrum))
      (is (= "test-001" (:id spectrum)))
      (is (= {:sample-type "test"} (:metadata spectrum)))
      (is (spectrum/dtype-spectrum? spectrum))))

  (testing "Spectrum creation with custom options"
    (let [[mz-vals intensities] (create-test-data 50)
          spectrum (core/create-spectrum "test-002" mz-vals intensities {}
                                         :validate? true
                                         :optimize? true
                                         :mz-dtype :float64
                                         :intensity-dtype :float32)]
      (is (spectrum/valid-spectrum? spectrum))
      (is (spectrum/dtype-spectrum? spectrum))))

  (testing "Legacy spectrum creation"
    (let [[mz-vals intensities] (create-test-data 50)
          spectrum (core/create-spectrum "test-003" mz-vals intensities {}
                                         :optimize? false)]
      (is (spectrum/valid-spectrum? spectrum))
      (is (spectrum/legacy-spectrum? spectrum))))

  (testing "Invalid spectrum creation"
    (is (thrown? Exception
                 (core/create-spectrum "test-004" [] [1 2 3] {})))
    (is (thrown? Exception
                 (core/create-spectrum "test-005" [1 2] [1 2 3] {})))))

(deftest test-spectrum-validation
  "Test spectrum validation functions."
  (testing "Valid spectrum detection"
    (let [[mz-vals intensities] (create-test-data 100)
          spectrum (core/create-spectrum "test" mz-vals intensities {})]
      (is (core/valid-spectrum? spectrum))))

  (testing "Invalid spectrum detection"
    (let [invalid-spectrum {:id "bad" :mz-values [] :intensities [1 2 3]}]
      (is (not (core/valid-spectrum? invalid-spectrum))))))

(deftest test-spectrum-info
  "Test spectrum information extraction."
  (testing "Comprehensive spectrum information"
    (let [[mz-vals intensities] (create-test-data 200)
          spectrum (core/create-spectrum "info-test" mz-vals intensities
                                         {:sample-type "E. coli"
                                          :date "2025-08-09"})
          info (core/spectrum-info spectrum)]
      (is (= "info-test" (:id info)))
      (is (= {:sample-type "E. coli" :date "2025-08-09"} (:metadata info)))
      (is (= 200 (:length info)))
      (is (vector? (:mz-range info)))
      (is (= 2 (count (:mz-range info))))
      (is (map? (:intensity-stats info)))
      (is (number? (:total-ion-current info)))
      (is (boolean? (:is-empty? info)))
      (is (boolean? (:dtype-optimized? info)))
      (is (not (:is-empty? info)))
      (is (:dtype-optimized? info))))

  (testing "Empty spectrum information"
    (let [empty-spectrum (spectrum/create-spectrum "empty" [] [] {})]
      (is (:is-empty? (core/spectrum-info empty-spectrum))))))

;; =============================================================================
;; Preprocessing Function Tests
;; =============================================================================

(deftest test-clinical-preprocess
  "Test clinical preprocessing pipeline."
  (testing "Default clinical preprocessing"
    (let [[mz-vals intensities] (create-test-data 1000)
          spectrum (core/create-spectrum "clinical-001" mz-vals intensities {})
          processed (core/clinical-preprocess spectrum)]
      (is (spectrum/valid-spectrum? processed))
      ;; Should be trimmed to clinical range
      (let [mz-range (:mz-range (core/spectrum-info processed))]
        (is (>= (first mz-range) 2000))
        (is (<= (second mz-range) 20000)))
      ;; Should be TIC normalized (sum ~= 1.0)
      (is (< (Math/abs (- (dfn/sum (:intensities processed)) 1.0)) 0.01))))

  (testing "Custom clinical preprocessing parameters"
    (let [[mz-vals intensities] (create-test-data 500)
          spectrum (core/create-spectrum "clinical-002" mz-vals intensities {})
          processed (core/clinical-preprocess spectrum
                                              :trim-range [1200 1800]
                                              :savgol-window 15
                                              :snip-iterations 50
                                              :tic-normalize? false)]
      (is (spectrum/valid-spectrum? processed))
      ;; Should be trimmed to custom range
      (let [mz-range (:mz-range (core/spectrum-info processed))]
        (is (>= (first mz-range) 1200))
        (is (<= (second mz-range) 1800)))
      ;; Should NOT be TIC normalized
      (is (> (dfn/sum (:intensities processed)) 10.0)))))

(deftest test-basic-preprocess
  "Test basic preprocessing pipeline."
  (testing "Default basic preprocessing"
    (let [[mz-vals intensities] (create-test-data 200)
          spectrum (core/create-spectrum "basic-001" mz-vals intensities {})
          processed (core/basic-preprocess spectrum)]
      (is (spectrum/valid-spectrum? processed))
      ;; Should be TIC normalized
      (is (< (Math/abs (- (dfn/sum (:intensities processed)) 1.0)) 0.01))))

  (testing "Custom basic preprocessing options"
    (let [[mz-vals intensities] (create-test-data 100)
          spectrum (core/create-spectrum "basic-002" mz-vals intensities {})
          processed (core/basic-preprocess spectrum
                                           :smooth? false
                                           :baseline-correct? true
                                           :baseline-method :snip
                                           :normalize? false)]
      (is (spectrum/valid-spectrum? processed))
      ;; Should NOT be normalized
      (is (> (dfn/sum (:intensities processed)) 10.0))))

  (testing "All preprocessing disabled"
    (let [[mz-vals intensities] (create-test-data 100)
          spectrum (core/create-spectrum "basic-003" mz-vals intensities {})
          processed (core/basic-preprocess spectrum
                                           :smooth? false
                                           :baseline-correct? false
                                           :normalize? false)]
      (is (spectrum/valid-spectrum? processed))
      ;; Should be essentially unchanged
      (is (= (spectrum/spectrum-length spectrum)
             (spectrum/spectrum-length processed))))))

;; =============================================================================
;; Peak Detection Tests
;; =============================================================================

(deftest test-detect-peaks
  "Test peak detection functionality."
  (testing "Basic peak detection with preprocessing"
    (let [[mz-vals intensities] (create-test-data 300 :with-peaks? true)
          spectrum (core/create-spectrum "peaks-001" mz-vals intensities {})
          peaks (core/detect-peaks spectrum)]
      (is (vector? peaks))
      (is (pos? (count peaks)))
      ;; Each peak should have required keys
      (doseq [peak peaks]
        (is (contains? peak :mz))
        (is (contains? peak :intensity))
        (is (contains? peak :snr)))))

  (testing "Peak detection without preprocessing"
    (let [[mz-vals intensities] (create-test-data 200 :with-peaks? true)
          spectrum (core/create-spectrum "peaks-002" mz-vals intensities {})
          peaks (core/detect-peaks spectrum :preprocess? false)]
      (is (vector? peaks))
      (is (pos? (count peaks)))))

  (testing "Custom peak detection parameters"
    (let [[mz-vals intensities] (create-test-data 200 :with-peaks? true)
          spectrum (core/create-spectrum "peaks-003" mz-vals intensities {})
          strict-peaks (core/detect-peaks spectrum
                                          :snr-threshold 5.0
                                          :window-size 5)
          lenient-peaks (core/detect-peaks spectrum
                                           :snr-threshold 1.0
                                           :window-size 15)]
      (is (<= (count strict-peaks) (count lenient-peaks)))))

  (testing "Peak detection on flat spectrum"
    (let [flat-intensities (repeat 100 50.0)
          spectrum (core/create-spectrum "flat" (range 1000 1100) flat-intensities {})
          peaks (core/detect-peaks spectrum)]
      ;; Should find few or no peaks in flat spectrum
      (is (<= (count peaks) 5)))))

(deftest test-peak-statistics
  "Test peak statistics calculation."
  (testing "Peak statistics calculation"
    (let [[mz-vals intensities] (create-test-data 300 :with-peaks? true)
          spectrum (core/create-spectrum "stats-001" mz-vals intensities {})
          peaks (core/detect-peaks spectrum)
          stats (core/peak-statistics peaks)]
      (is (map? stats))
      (is (contains? stats :count))
      (is (= (count peaks) (:count stats))))))

(deftest test-compare-spectra
  "Test spectrum comparison (placeholder implementation)."
  (testing "Basic spectrum comparison"
    (let [[mz-vals1 intensities1] (create-test-data 100)
          [mz-vals2 intensities2] (create-test-data 100)
          spectrum1 (core/create-spectrum "comp-1" mz-vals1 intensities1 {})
          spectrum2 (core/create-spectrum "comp-2" mz-vals2 intensities2 {})
          result (core/compare-spectra spectrum1 spectrum2)]
      (is (map? result))
      (is (= "placeholder" (:status result)))
      (is (contains? result :spectrum1-info))
      (is (contains? result :spectrum2-info)))))

;; =============================================================================
;; Batch Processing Tests
;; =============================================================================

(deftest test-batch-processing
  "Test batch processing functions."
  (testing "Basic batch preprocessing"
    (let [test-spectra (repeatedly 5
                                   (fn []
                                     (let [[mz-vals intensities] (create-test-data 100)]
                                       (core/create-spectrum (str "batch-" (rand-int 1000))
                                                             mz-vals intensities {}))))
          processed (core/batch-preprocess test-spectra)]
      (is (= 5 (count processed)))
      (is (every? spectrum/valid-spectrum? processed))))

  (testing "Batch clinical preprocessing"
    (let [test-spectra (repeatedly 3
                                   (fn []
                                     (let [[mz-vals intensities] (create-test-data 200)]
                                       (core/create-spectrum (str "clinical-batch-" (rand-int 1000))
                                                             mz-vals intensities {}))))
          processed (core/batch-preprocess test-spectra :method :clinical)]
      (is (= 3 (count processed)))
      (is (every? spectrum/valid-spectrum? processed))))

  (testing "Batch peak detection"
    (let [test-spectra (repeatedly 3
                                   (fn []
                                     (let [[mz-vals intensities] (create-test-data 150 :with-peaks? true)]
                                       (core/create-spectrum (str "peak-batch-" (rand-int 1000))
                                                             mz-vals intensities {}))))
          peak-results (core/batch-detect-peaks test-spectra)]
      (is (= 3 (count peak-results)))
      (is (every? vector? peak-results))
      (is (every? #(pos? (count %)) peak-results)))))

;; =============================================================================
;; Data Export and Conversion Tests  
;; =============================================================================

(deftest test-data-conversion
  "Test data conversion functions."
  (testing "Spectra to dataset conversion"
    (let [test-spectra (repeatedly 3
                                   (fn []
                                     (let [[mz-vals intensities] (create-test-data 50)]
                                       (core/create-spectrum (str "dataset-" (rand-int 1000))
                                                             mz-vals intensities
                                                             {:sample-type "test"}))))
          dataset (core/spectra->dataset test-spectra)]
      (is (ds/dataset? dataset))
      (is (= 3 (ds/row-count dataset)))))

  (testing "Spectra to dataset with peaks"
    (let [test-spectra (repeatedly 2
                                   (fn []
                                     (let [[mz-vals intensities] (create-test-data 100 :with-peaks? true)]
                                       (core/create-spectrum (str "peaks-dataset-" (rand-int 1000))
                                                             mz-vals intensities {}))))
          dataset (core/spectra->dataset test-spectra :include-peaks? true)]
      (is (ds/dataset? dataset))
      (is (= 2 (ds/row-count dataset)))
      ;; Should have peak-related columns
      (is (some #(re-find #"peak" (str %)) (ds/column-names dataset))))))

;; =============================================================================
;; Utility Function Tests
;; =============================================================================

(deftest test-version-info
  "Test system information function."
  (testing "Version information structure"
    (let [info (core/version-info)]
      (is (map? info))
      (is (contains? info :maldicloj-version))
      (is (contains? info :clojure-version))
      (is (contains? info :java-version))
      (is (contains? info :available-memory-mb))
      (is (contains? info :processors))
      (is (string? (:maldicloj-version info)))
      (is (string? (:clojure-version info)))
      (is (number? (:available-memory-mb info)))
      (is (number? (:processors info))))))

(deftest test-benchmark-processing
  "Test performance benchmarking."
  (testing "Small benchmark test"
    (let [results (core/benchmark-processing 3 50)]
      (is (map? results))
      (is (= 3 (:n-spectra results)))
      (is (= 50 (:n-points results)))
      (is (contains? results :clinical-preprocessing))
      (is (contains? results :peak-detection))
      (is (map? (:clinical-preprocessing results)))
      (is (map? (:peak-detection results)))
      (is (number? (get-in results [:clinical-preprocessing :total-time-ms])))
      (is (vector? (:peak-counts results)))
      (is (= 3 (count (:peak-counts results)))))))

(deftest test-help-function
  "Test help display function."
  (testing "Help function executes without error"
    (is (nil? (core/help))))) ; help returns nil after printing

;; =============================================================================
;; Error Handling and Edge Cases
;; =============================================================================

(deftest test-error-handling
  "Test error handling in core functions."
  (testing "Invalid spectrum input handling"
    (let [invalid-spectrum {:not-a-spectrum true}]
      (is (thrown? Exception (core/clinical-preprocess invalid-spectrum)))
      (is (thrown? Exception (core/detect-peaks invalid-spectrum)))))

  (testing "Empty batch processing"
    (let [empty-results (core/batch-preprocess [])]
      (is (empty? empty-results))))

  (testing "Export format error handling"
    (is (thrown? Exception
                 (core/export-results {:data "test"} :invalid-format "test.txt")))))

;; =============================================================================
;; Integration Tests
;; =============================================================================

(deftest test-complete-workflow
  "Test complete analysis workflow."
  (testing "End-to-end spectrum analysis"
    (let [[mz-vals intensities] (create-test-data 300 :with-peaks? true)
          ;; Step 1: Create spectrum
          spectrum (core/create-spectrum "workflow-test" mz-vals intensities
                                         {:sample-type "E. coli"
                                          :operator "test-user"})

          ;; Step 2: Get spectrum info
          info (core/spectrum-info spectrum)

          ;; Step 3: Clinical preprocessing
          processed (core/clinical-preprocess spectrum)

          ;; Step 4: Peak detection
          peaks (core/detect-peaks processed :preprocess? false)

          ;; Step 5: Peak statistics
          peak-stats (core/peak-statistics peaks)]

      ;; Validate each step
      (is (spectrum/valid-spectrum? spectrum))
      (is (map? info))
      (is (:dtype-optimized? info))
      (is (spectrum/valid-spectrum? processed))
      (is (vector? peaks))
      (is (pos? (count peaks)))
      (is (map? peak-stats))
      (is (= (count peaks) (:count peak-stats))))))

;; =============================================================================
;; Property-Based Tests
;; =============================================================================

(defspec test-spectrum-creation-properties 25
  (prop/for-all [n-points (gen/choose 10 200)
                 id (gen/not-empty gen/string-alphanumeric)]
                (let [[mz-vals intensities] (create-test-data n-points)
                      spectrum (core/create-spectrum id mz-vals intensities {})]
                  (and (spectrum/valid-spectrum? spectrum)
                       (= id (:id spectrum))
                       (= n-points (spectrum/spectrum-length spectrum))))))

(defspec test-preprocessing-preserves-validity 20
  (prop/for-all [n-points (gen/choose 50 150)]
                (let [[mz-vals intensities] (create-test-data n-points)
                      spectrum (core/create-spectrum "prop-test" mz-vals intensities {})
                      basic-processed (core/basic-preprocess spectrum)
                      clinical-processed (core/clinical-preprocess spectrum)]
                  (and (spectrum/valid-spectrum? basic-processed)
                       (spectrum/valid-spectrum? clinical-processed)))))

(defspec test-peak-detection-consistency 15
  (prop/for-all [n-points (gen/choose 100 200)
                 snr-threshold (gen/choose 1.0 5.0)]
                (let [[mz-vals intensities] (create-test-data n-points :with-peaks? true)
                      spectrum (core/create-spectrum "consistency-test" mz-vals intensities {})
                      peaks1 (core/detect-peaks spectrum :snr-threshold snr-threshold)
                      peaks2 (core/detect-peaks spectrum :snr-threshold snr-threshold)]
      ;; Same parameters should give same results
                  (= (count peaks1) (count peaks2)))))

;; =============================================================================
;; Test Runner
;; =============================================================================

(defn run-core-tests
  "Run all core API tests and return results."
  []
  (run-tests 'maldi-clj.core-test))
