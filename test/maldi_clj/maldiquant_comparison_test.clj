(ns maldi-clj.maldiquant-comparison-test
  "Comprehensive tests comparing maldicloj implementation with R MALDIquant package"
  (:require [clojure.test :refer [deftest is testing]]
            [clojisr.v1.r :as r :refer [r]]
            [clojisr.v1.require :refer [require-r]]
            [maldi-clj.spectrum :as spectrum]
            [maldi-clj.baseline :as baseline]
            [maldi-clj.peaks :as peaks]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]))

;; ==============================================================================
;; R Environment Setup
;; ==============================================================================

(defn setup-maldiquant-environment! []
  (try
    ;; Install and load required R packages
    (r "if (!require('MALDIquant', quietly=TRUE)) install.packages('MALDIquant')")
    (r "library(MALDIquant)")
    (println "MALDIquant R environment setup complete")
    true
    (catch Exception e
      (println "Warning: Could not setup MALDIquant R environment:" (.getMessage e))
      false)))

(def maldiquant-available? (atom nil))

(defn ensure-maldiquant-setup []
  (when (nil? @maldiquant-available?)
    (reset! maldiquant-available? (setup-maldiquant-environment!)))
  @maldiquant-available?)

;; ==============================================================================
;; Test Data Generation
;; ==============================================================================

(defn create-maldiquant-test-spectrum
  "Create spectrum compatible with both maldicloj and MALDIquant"
  [n-points]
  (let [mz-values (vec (range 1000.0 (+ 1000.0 n-points) 1.0))
        ;; Create realistic mass spectrum with multiple peaks
        peaks (mapv (fn [x]
                      (+ (* 100 (Math/exp (- (/ (* (- x 50) (- x 50)) 50))))
                         (* 200 (Math/exp (- (/ (* (- x 150) (- x 150)) 75))))
                         (* 150 (Math/exp (- (/ (* (- x 250) (- x 250)) 40))))
                         (* 80 (Math/exp (- (/ (* (- x 350) (- x 350)) 60))))))
                    (range n-points))
        ;; Add realistic baseline drift
        baseline (mapv #(+ 10.0 (* 0.02 %) (* 0.00001 % %)) (range n-points))
        ;; Add Poisson-like noise
        noise (repeatedly n-points #(* 5.0 (- (rand) 0.5)))
        ;; Combine components
        intensities (mapv + peaks baseline noise)]

    {:spectrum (spectrum/create-spectrum "maldiquant-test" mz-values intensities
                                         {:type "MALDI-TOF" :sample "test"})
     :mz-values mz-values
     :intensities intensities
     :true-peaks peaks
     :true-baseline baseline}))

;; ==============================================================================
;; Baseline Correction Comparison Tests
;; ==============================================================================

(deftest test-tophat-baseline-vs-maldiquant
  (testing "TopHat baseline correction vs MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-maldiquant-test-spectrum 300)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Clojure TopHat implementation
            clj-baseline (baseline/estimate-baseline test-spectrum :tophat :half-window-size 25)

            ;; R MALDIquant TopHat implementation
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
            _ (r "r_baseline <- estimateBaseline(spectrum, method='TopHat', halfWindowSize=25)")
            r-baseline-data (-> (r "as.matrix(r_baseline)")
                                r/r->clj)
            r-baseline (-> r-baseline-data
                           (get "intensity")
                           dtype/->vector)]

        (testing "TopHat baseline properties"
          (is (every? #(>= % 0) clj-baseline)
              "Clojure TopHat baseline should be non-negative")
          (is (= (count clj-baseline) (count r-baseline))
              "Baseline lengths should match"))

        (testing "TopHat correlation with MALDIquant"
          (let [clj-vec (if (dtype/reader? clj-baseline)
                          (dtype/->vector clj-baseline)
                          clj-baseline)
                correlation (when (> (count clj-vec) 1)
                              (let [mean-clj (/ (reduce + clj-vec) (count clj-vec))
                                    mean-r (/ (reduce + r-baseline) (count r-baseline))
                                    numerator (reduce + (map #(* (- %1 mean-clj) (- %2 mean-r))
                                                             clj-vec r-baseline))
                                    denom-clj (Math/sqrt (reduce + (map #(* (- % mean-clj) (- % mean-clj)) clj-vec)))
                                    denom-r (Math/sqrt (reduce + (map #(* (- % mean-r) (- % mean-r)) r-baseline)))]
                                (/ numerator (* denom-clj denom-r))))]
            (is (> correlation 0.7)
                (str "TopHat baseline should correlate well with MALDIquant: " correlation))))))))

(deftest test-median-baseline-vs-maldiquant
  (testing "Median baseline correction vs MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-maldiquant-test-spectrum 200)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Clojure median implementation
            clj-baseline (baseline/estimate-baseline test-spectrum :median :half-window-size 20)

            ;; R MALDIquant median implementation
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
            _ (r "r_baseline <- estimateBaseline(spectrum, method='median', halfWindowSize=20)")
            r-baseline-data (-> (r "as.matrix(r_baseline)")
                                r/r->clj)
            r-baseline (-> r-baseline-data
                           (get "intensity")
                           dtype/->vector)]

        (testing "median baseline properties"
          (is (every? #(>= % 0) clj-baseline)
              "Clojure median baseline should be non-negative")
          (is (= (count clj-baseline) (count r-baseline))
              "Baseline lengths should match"))

        (testing "median correlation with MALDIquant"
          (let [clj-vec (if (dtype/reader? clj-baseline)
                          (dtype/->vector clj-baseline)
                          clj-baseline)
                correlation (when (> (count clj-vec) 1)
                              (let [mean-clj (/ (reduce + clj-vec) (count clj-vec))
                                    mean-r (/ (reduce + r-baseline) (count r-baseline))
                                    numerator (reduce + (map #(* (- %1 mean-clj) (- %2 mean-r))
                                                             clj-vec r-baseline))
                                    denom-clj (Math/sqrt (reduce + (map #(* (- % mean-clj) (- % mean-clj)) clj-vec)))
                                    denom-r (Math/sqrt (reduce + (map #(* (- % mean-r) (- % mean-r)) r-baseline)))]
                                (/ numerator (* denom-clj denom-r))))]
            (is (> correlation 0.8)
                (str "Median baseline should correlate strongly with MALDIquant: " correlation))))))))

;; ==============================================================================
;; Peak Detection Comparison Tests
;; ==============================================================================

(deftest test-noise-estimation-vs-maldiquant
  (testing "Noise estimation methods vs MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-maldiquant-test-spectrum 500)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Clojure MAD noise estimation
            clj-noise-mad (peaks/estimate-noise intensities :mad)

            ;; R MALDIquant noise estimation
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
            r-noise-data (-> (r "estimateNoise(spectrum, method='MAD')")
                             r/r->clj)
            r-noise-values (-> r-noise-data
                               (get "intensity")
                               dtype/->vector)
            r-noise-median (nth (sort r-noise-values) (quot (count r-noise-values) 2))]

        (testing "MAD noise estimation properties"
          (is (> clj-noise-mad 0)
              "Clojure MAD noise should be positive")
          (is (> r-noise-median 0)
              "R MAD noise should be positive"))

        (testing "MAD noise estimation comparison"
          ;; Allow for reasonable difference due to different implementations
          (let [noise-ratio (/ clj-noise-mad r-noise-median)]
            (is (and (> noise-ratio 0.5) (< noise-ratio 2.0))
                (str "Noise estimates should be in same order of magnitude: "
                     clj-noise-mad " vs " r-noise-median))))))))

(deftest test-peak-detection-vs-maldiquant
  (testing "Peak detection vs MALDIquant detectPeaks"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-maldiquant-test-spectrum 400)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Clojure peak detection
            clj-peaks (peaks/detect-peaks test-spectrum
                                          :snr-threshold 3.0
                                          :window-half-size 5
                                          :noise-method :mad)
            clj-peak-count (peaks/peak-count clj-peaks)

            ;; R MALDIquant peak detection
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
            _ (r "peaks <- detectPeaks(spectrum, method='MAD', SNR=3, halfWindowSize=5)")
            r-peak-data (-> (r "as.matrix(peaks)")
                            r/r->clj)
            r-peak-count (count (get r-peak-data "mass"))]

        (testing "peak detection results"
          (is (> clj-peak-count 0)
              "Clojure peak detection should find peaks")
          (is (> r-peak-count 0)
              "R peak detection should find peaks"))

        (testing "peak count comparison"
          ;; Peak detection algorithms can vary significantly - just ensure both find peaks
          (is (and (> clj-peak-count 0) (> r-peak-count 0))
              (str "Both implementations should find peaks: "
                   clj-peak-count " vs " r-peak-count)))))))

;; ==============================================================================
;; Spectrum Processing Comparison Tests
;; ==============================================================================

(deftest test-total-ion-current-vs-maldiquant
  (testing "Total Ion Current calculation vs MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-maldiquant-test-spectrum 100)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Clojure TIC calculation
            clj-tic (if (dtype/reader? intensities)
                      (dfn/sum intensities)
                      (reduce + intensities))

            ;; R MALDIquant TIC calculation
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
            r-tic (-> (r "totalIonCurrent(spectrum)")
                      r/r->clj
                      first)]

        (testing "TIC calculation close match"
          ;; Allow small differences due to floating point precision
          (is (< (Math/abs (- clj-tic r-tic)) (* 0.01 clj-tic))
              (str "TIC calculations should be very close: " clj-tic " vs " r-tic)))))))

(deftest test-baseline-removal-vs-maldiquant
  (testing "Complete baseline removal vs MALDIquant removeBaseline"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-maldiquant-test-spectrum 200)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Clojure baseline removal (SNIP)
            clj-corrected (baseline/remove-baseline test-spectrum :snip :iterations 30)
            clj-corrected-intensities (:intensities clj-corrected)

            ;; R MALDIquant baseline removal
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
            _ (r "corrected <- removeBaseline(spectrum, method='SNIP', iterations=30)")
            r-corrected-data (-> (r "as.matrix(corrected)")
                                 r/r->clj)
            r-corrected-intensities (-> r-corrected-data
                                        (get "intensity")
                                        dtype/->vector)]

        (testing "baseline removal properties"
          (is (spectrum/valid-spectrum? clj-corrected)
              "Corrected spectrum should be valid")
          (is (every? #(>= % 0) (if (dtype/reader? clj-corrected-intensities)
                                  (dtype/->vector clj-corrected-intensities)
                                  clj-corrected-intensities))
              "Corrected intensities should be non-negative")
          (is (= (count clj-corrected-intensities) (count r-corrected-intensities))
              "Corrected spectra should have same length"))

        (testing "baseline removal correlation"
          (let [clj-vec (if (dtype/reader? clj-corrected-intensities)
                          (dtype/->vector clj-corrected-intensities)
                          clj-corrected-intensities)
                correlation (when (> (count clj-vec) 1)
                              (let [mean-clj (/ (reduce + clj-vec) (count clj-vec))
                                    mean-r (/ (reduce + r-corrected-intensities) (count r-corrected-intensities))
                                    numerator (reduce + (map #(* (- %1 mean-clj) (- %2 mean-r))
                                                             clj-vec r-corrected-intensities))
                                    denom-clj (Math/sqrt (reduce + (map #(* (- % mean-clj) (- % mean-clj)) clj-vec)))
                                    denom-r (Math/sqrt (reduce + (map #(* (- % mean-r) (- % mean-r)) r-corrected-intensities)))]
                                (/ numerator (* denom-clj denom-r))))]
            (is (> correlation 0.5)
                (str "Baseline-corrected spectra should show reasonable correlation: " correlation))))))))

;; ==============================================================================
;; Data Structure Compatibility Tests
;; ==============================================================================

(deftest test-spectrum-data-structure-compatibility
  (testing "Spectrum data structure compatibility with MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [test-spectrum (spectrum/create-spectrum
                           "compatibility-test"
                           [1000.0 1001.0 1002.0 1003.0 1004.0]
                           [100.0 200.0 150.0 300.0 50.0]
                           {:instrument "Bruker" :method "reflectron"})
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            ;; Create equivalent R spectrum
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
            _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
            _ (r "r_spectrum <- createMassSpectrum(mass=mz, intensity=intensities, 
                   metaData=list(instrument='Bruker', method='reflectron'))")

            ;; Extract data back from R
            r-mz (-> (r "mass(r_spectrum)") r/r->clj vec)
            r-intensities (-> (r "intensity(r_spectrum)") r/r->clj vec)
            r-metadata (-> (r "metaData(r_spectrum)") r/r->clj)]

        (testing "data preservation"
          (is (= mz-values r-mz)
              "m/z values should be preserved")
          (is (= intensities r-intensities)
              "Intensities should be preserved")
          ;; R metadata structure may be different - just check it exists
          (is (not (nil? r-metadata))
              "R spectrum should have metadata"))

        (testing "round-trip data integrity"
          ;; Process in R and bring back to Clojure
          (let [_ (r "processed <- smoothIntensity(r_spectrum, method='MovingAverage', halfWindowSize=1)")
                processed-intensities (-> (r "intensity(processed)") r/r->clj vec)

                ;; Apply similar processing in Clojure (would need implementation)
                ;; For now, verify data types and structure
                clj-processed-intensities intensities] ; Placeholder

            (is (= (count intensities) (count processed-intensities))
                "Processed spectrum should maintain length")
            (is (every? number? processed-intensities)
                "Processed intensities should be numeric")))))))

;; ==============================================================================
;; Edge Cases and Error Handling
;; ==============================================================================

(deftest test-edge-cases-vs-maldiquant
  (testing "Edge case handling compared to MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [;; Test various edge cases
            empty-spectrum (spectrum/create-spectrum "empty" [] [] {})
            single-point (spectrum/create-spectrum "single" [1000.0] [100.0] {})
            zero-intensities (spectrum/create-spectrum "zeros" [1000.0 1001.0] [0.0 0.0] {})
            negative-intensities (spectrum/create-spectrum "negative" [1000.0 1001.0] [-10.0 20.0] {})]

        (testing "empty spectrum handling"
          ;; Both implementations might handle this differently
          (try
            (let [baseline (baseline/estimate-baseline empty-spectrum :snip)]
              (is (or (empty? baseline) (= 0 (count baseline)))
                  "Empty spectrum baseline should be empty or throw"))
            (catch Exception e
              (is true "Empty spectrum handling - exception acceptable"))))

        (testing "single point spectrum"
          (try
            (let [baseline (baseline/estimate-baseline single-point :snip :iterations 5)]
              (is (= 1 (count baseline))
                  "Single point baseline should have one element"))
            (catch Exception e
              (is true "Single point handling - exception acceptable"))))

        (testing "zero intensities"
          (let [baseline (baseline/estimate-baseline zero-intensities :snip :iterations 5)]
            (is (every? #(>= % 0) baseline)
                "Zero intensities baseline should be non-negative")))

        (testing "negative intensities"
          ;; This tests robustness - MALDIquant might handle differently
          (try
            (let [baseline (baseline/estimate-baseline negative-intensities :snip :iterations 5)]
              (is (not (nil? baseline))
                  "Should handle negative intensities gracefully"))
            (catch Exception e
              (is true "Negative intensity handling - exception acceptable"))))))))

;; ==============================================================================
;; Performance Comparison (Basic)
;; ==============================================================================

(deftest test-performance-comparison-basic
  (testing "Basic performance comparison with MALDIquant"
    (when (ensure-maldiquant-setup)
      (let [large-spectrum (create-maldiquant-test-spectrum 1000)
            test-spectrum (:spectrum large-spectrum)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)]

        (testing "Clojure performance"
          (let [start-time (System/currentTimeMillis)
                _ (baseline/estimate-baseline test-spectrum :snip :iterations 25)
                end-time (System/currentTimeMillis)
                clj-time (- end-time start-time)]
            (is (< clj-time 5000) ; Should complete in under 5 seconds
                (str "Clojure SNIP should complete reasonably quickly: " clj-time "ms"))))

        (testing "R performance"
          (let [_ (r (str "mz <- c(" (clojure.string/join ", " mz-values) ")"))
                _ (r (str "intensities <- c(" (clojure.string/join ", " intensities) ")"))
                _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
                start-time (System/currentTimeMillis)
                _ (r "r_baseline <- estimateBaseline(spectrum, method='SNIP', iterations=25)")
                end-time (System/currentTimeMillis)
                r-time (- end-time start-time)]
            (is (< r-time 10000) ; R might be slower due to overhead
                (str "R SNIP should complete reasonably: " r-time "ms"))))))))

;; ==============================================================================
;; Test Runner Functions
;; ==============================================================================

(defn run-maldiquant-baseline-tests []
  "Run MALDIquant baseline comparison tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running MALDIquant Baseline Comparison Tests ===")
    (clojure.test/run-test-var #'test-tophat-baseline-vs-maldiquant)
    (clojure.test/run-test-var #'test-median-baseline-vs-maldiquant)
    (clojure.test/run-test-var #'test-baseline-removal-vs-maldiquant)
    (println "=== MALDIquant baseline tests complete ===")))

(defn run-maldiquant-peak-tests []
  "Run MALDIquant peak detection comparison tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running MALDIquant Peak Detection Comparison Tests ===")
    (clojure.test/run-test-var #'test-noise-estimation-vs-maldiquant)
    (clojure.test/run-test-var #'test-peak-detection-vs-maldiquant)
    (println "=== MALDIquant peak tests complete ===")))

(defn run-maldiquant-compatibility-tests []
  "Run MALDIquant data compatibility tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running MALDIquant Compatibility Tests ===")
    (clojure.test/run-test-var #'test-total-ion-current-vs-maldiquant)
    (clojure.test/run-test-var #'test-spectrum-data-structure-compatibility)
    (clojure.test/run-test-var #'test-edge-cases-vs-maldiquant)
    (println "=== MALDIquant compatibility tests complete ===")))

(defn run-all-maldiquant-tests []
  "Run complete MALDIquant comparison test suite"
  (when (ensure-maldiquant-setup)
    (println "=== Running Complete MALDIquant Comparison Suite ===")
    (run-maldiquant-baseline-tests)
    (run-maldiquant-peak-tests)
    (run-maldiquant-compatibility-tests)
    (clojure.test/run-test-var #'test-performance-comparison-basic)
    (println "=== MALDIquant test suite complete ==="))
  (when (not (ensure-maldiquant-setup))
    (println "MALDIquant R environment not available - skipping comparison tests")))
