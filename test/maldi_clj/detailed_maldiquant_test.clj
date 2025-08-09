(ns maldi-clj.detailed-maldiquant-test
  "Comprehensive detailed tests comparing maldicloj with R MALDIquant package"
  (:require [clojure.test :refer [deftest is testing]]
            [clojisr.v1.r :as r :refer [r]]
            [maldi-clj.spectrum :as spectrum]
            [maldi-clj.baseline :as baseline]
            [maldi-clj.peaks :as peaks]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [clojure.string :as str]))

;; ==============================================================================
;; Enhanced Test Data Generation
;; ==============================================================================

(defn create-realistic-maldi-spectrum
  "Create realistic MALDI-TOF spectrum with controlled characteristics"
  [n-points & {:keys [peak-positions peak-intensities baseline-drift noise-level
                      mz-range instrument-resolution]
               :or {peak-positions [1200 1800 2400 3000 3600]
                    peak-intensities [1000 800 1200 600 400]
                    baseline-drift 0.02
                    noise-level 10.0
                    mz-range [1000 4000]
                    instrument-resolution 1000}}]
  (let [mz-start (first mz-range)
        mz-end (second mz-range)
        mz-values (vec (range mz-start mz-end (/ (- mz-end mz-start) n-points)))

        ;; Generate realistic peak shapes
        peak-signals (map (fn [pos intensity]
                            (mapv (fn [mz]
                                    (let [sigma (/ pos instrument-resolution)
                                          normalized-dist (/ (- mz pos) sigma)
                                          gaussian (* intensity (Math/exp (- (/ (* normalized-dist normalized-dist) 2))))
                                          asymmetry-factor (if (> mz pos) 0.9 1.0)]
                                      (* gaussian asymmetry-factor)))
                                  mz-values))
                          peak-positions peak-intensities)

        ;; Combine all peaks
        combined-peaks (reduce (fn [acc peak-signal]
                                 (mapv + acc peak-signal))
                               (vec (repeat n-points 0.0))
                               peak-signals)

        ;; Add realistic baseline drift
        baseline (mapv #(+ 5.0 (* baseline-drift %) (* 0.000005 % %)) mz-values)

        ;; Add realistic noise
        noise (repeatedly n-points #(+ (* noise-level (- (rand) 0.5))
                                       (* 2.0 (Math/sqrt (max 1 (rand-int 50))))))

        ;; Combine all components
        intensities (mapv #(max 0.1 (+ %1 %2 %3)) combined-peaks baseline noise)]

    {:spectrum (spectrum/create-spectrum "realistic-maldi" mz-values intensities
                                         {:instrument "MALDI-TOF" :resolution instrument-resolution})
     :mz-values mz-values
     :intensities intensities
     :true-peaks combined-peaks
     :true-baseline baseline
     :peak-positions peak-positions
     :peak-intensities peak-intensities
     :noise noise}))

;; ==============================================================================
;; Setup and Utilities
;; ==============================================================================

(defn ensure-maldiquant-setup []
  (try
    (r "library(MALDIquant)")
    true
    (catch Exception e
      (println "MALDIquant not available:" (.getMessage e))
      false)))

(defn calculate-correlation
  "Calculate Pearson correlation between two sequences"
  [seq1 seq2]
  (when (and (> (count seq1) 1) (= (count seq1) (count seq2)))
    (let [n (count seq1)
          mean1 (/ (reduce + seq1) n)
          mean2 (/ (reduce + seq2) n)
          numerator (reduce + (map #(* (- %1 mean1) (- %2 mean2)) seq1 seq2))
          denom1 (Math/sqrt (reduce + (map #(* (- % mean1) (- % mean1)) seq1)))
          denom2 (Math/sqrt (reduce + (map #(* (- % mean2) (- % mean2)) seq2)))]
      (if (and (> denom1 0) (> denom2 0))
        (/ numerator (* denom1 denom2))
        0.0))))

(defn to-r-vector [data]
  "Convert Clojure sequence to R vector string"
  (str "c(" (str/join ", " data) ")"))

(defn statistical-summary [data]
  "Calculate statistical summary of data"
  (let [n (count data)
        mean (/ (reduce + data) n)
        variance (/ (reduce + (map #(* (- % mean) (- % mean)) data)) n)
        std-dev (Math/sqrt variance)]
    {:mean mean :std-dev std-dev :variance variance :count n}))

;; ==============================================================================
;; Detailed Baseline Correction Tests
;; ==============================================================================

(deftest test-snip-parameter-sensitivity-detailed
  (testing "Comprehensive SNIP parameter sensitivity analysis"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-realistic-maldi-spectrum 400)
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            iteration-counts [10 25 50 100 150]

            ;; Setup R spectrum once
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- " (to-r-vector mz-values)))
            _ (r (str "intensities <- " (to-r-vector intensities)))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")]

        (testing "iteration count sensitivity"
          (doseq [iter-count iteration-counts]
            (let [clj-result (baseline/estimate-baseline test-spectrum :snip :iterations iter-count)
                  _ (r (str "r_baseline <- estimateBaseline(spectrum, method='SNIP', iterations=" iter-count ")"))
                  r-result (-> (r "as.vector(r_baseline)[-(1:length(mz))]") r/r->clj vec)

                  clj-vec (if (dtype/reader? clj-result) (dtype/->vector clj-result) clj-result)
                  correlation (calculate-correlation clj-vec r-result)]

              (is (> correlation 0.5)
                  (str "SNIP correlation should be reasonable for " iter-count " iterations: " correlation))

              (is (every? #(>= % 0) clj-vec)
                  (str "SNIP baseline should be non-negative for " iter-count " iterations")))))

        (testing "parameter effect on baseline quality"
          (let [low-iter-baseline (baseline/estimate-baseline test-spectrum :snip :iterations 10)
                high-iter-baseline (baseline/estimate-baseline test-spectrum :snip :iterations 100)

                low-iter-vec (if (dtype/reader? low-iter-baseline) (dtype/->vector low-iter-baseline) low-iter-baseline)
                high-iter-vec (if (dtype/reader? high-iter-baseline) (dtype/->vector high-iter-baseline) high-iter-baseline)

                low-stats (statistical-summary low-iter-vec)
                high-stats (statistical-summary high-iter-vec)]

            ;; High iterations should generally produce smoother baseline
            (is (<= (:std-dev high-stats) (* 1.2 (:std-dev low-stats)))
                "More iterations should produce similar or smoother baseline")))))))

(deftest test-baseline-methods-comprehensive-detailed
  (testing "Comprehensive comparison of baseline methods with detailed analysis"
    (when (ensure-maldiquant-setup)
      (let [test-scenarios [{:name "low-noise" :data (create-realistic-maldi-spectrum 300 :noise-level 5.0)}
                            {:name "medium-noise" :data (create-realistic-maldi-spectrum 300 :noise-level 15.0)}
                            {:name "high-noise" :data (create-realistic-maldi-spectrum 300 :noise-level 30.0)}]

            methods [:snip :tophat :median]
            r-method-names ["SNIP" "TopHat" "median"]]

        (doseq [{:keys [name data]} test-scenarios]
          (let [test-spectrum (:spectrum data)
                mz-values (:mz-values test-spectrum)
                intensities (:intensities test-spectrum)]

            (testing (str "scenario: " name)
              ;; Setup R spectrum
              (r "library(MALDIquant)")
              (r (str "mz <- " (to-r-vector mz-values)))
              (r (str "intensities <- " (to-r-vector intensities)))
              (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")

              (doseq [[method r-method] (map vector methods r-method-names)]
                (testing (str "method: " method)
                  (let [clj-baseline (baseline/estimate-baseline test-spectrum method
                                                                 :iterations 50 :half-window-size 25)

                        ;; R implementation
                        _ (case method
                            :snip (r "r_baseline <- estimateBaseline(spectrum, method='SNIP', iterations=50)")
                            :tophat (r "r_baseline <- estimateBaseline(spectrum, method='TopHat', halfWindowSize=25)")
                            :median (r "r_baseline <- estimateBaseline(spectrum, method='median', halfWindowSize=25)"))
                        r-baseline (-> (r "as.vector(r_baseline)[-(1:length(mz))]") r/r->clj vec)

                        clj-vec (if (dtype/reader? clj-baseline) (dtype/->vector clj-baseline) clj-baseline)
                        correlation (calculate-correlation clj-vec r-baseline)

                        clj-stats (statistical-summary clj-vec)
                        r-stats (statistical-summary r-baseline)]

                    (is (> correlation 0.4)
                        (str method " on " name " correlation: " correlation))

                    (is (every? #(>= % 0) clj-vec)
                        (str method " on " name " should be non-negative"))

                    (is (< (:std-dev clj-stats) (:mean clj-stats))
                        (str method " on " name " should be smoother than original"))

                    ;; Cross-implementation mean similarity
                    (is (< (Math/abs (- (:mean clj-stats) (:mean r-stats)))
                           (* 0.4 (:mean clj-stats)))
                        (str method " on " name " should have similar mean to R"))))))))))))

(deftest test-peak-detection-parameter-analysis
  (testing "Detailed peak detection parameter sensitivity"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-realistic-maldi-spectrum 400
                                                       :peak-positions [1200 1600 2000 2400 2800]
                                                       :peak-intensities [800 600 1000 400 500])
            test-spectrum (:spectrum test-data)
            expected-peaks (:peak-positions test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)

            snr-thresholds [1.5 2.0 3.0 5.0 8.0]

            ;; Setup R spectrum
            _ (r "library(MALDIquant)")
            _ (r (str "mz <- " (to-r-vector mz-values)))
            _ (r (str "intensities <- " (to-r-vector intensities)))
            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")]

        (testing "SNR threshold effects on peak count"
          (let [results (map (fn [snr-thresh]
                               (let [clj-peaks (peaks/detect-peaks test-spectrum
                                                                   :snr-threshold snr-thresh
                                                                   :window-half-size 8
                                                                   :noise-method :mad)
                                     clj-count (peaks/peak-count clj-peaks)]
                                 {:snr snr-thresh :count clj-count}))
                             snr-thresholds)]

            ;; Higher SNR should generally find fewer peaks
            (doseq [[low high] (partition 2 1 results)]
              (is (>= (:count low) (:count high))
                  (str "Higher SNR should find <= peaks: SNR " (:snr low)
                       " found " (:count low) ", SNR " (:snr high) " found " (:count high))))

            ;; Should find some peaks at reasonable thresholds
            (let [moderate-results (filter #(and (>= (:snr %) 2.0) (<= (:snr %) 5.0)) results)]
              (is (every? #(> (:count %) 0) moderate-results)
                  "Should find peaks at moderate SNR thresholds"))))

        (testing "peak localization quality"
          (let [clj-peaks (peaks/detect-peaks test-spectrum
                                              :snr-threshold 2.5
                                              :window-half-size 8
                                              :noise-method :mad)
                detected-mz (peaks/get-peak-mz-values clj-peaks)

                ;; Find closest match for each expected peak
                matches (map (fn [expected]
                               (if (seq detected-mz)
                                 (let [closest (apply min-key #(Math/abs (- % expected)) detected-mz)
                                       distance (Math/abs (- closest expected))]
                                   {:expected expected :detected closest :distance distance :found true})
                                 {:expected expected :found false}))
                             expected-peaks)

                close-matches (filter #(and (:found %) (< (:distance %) 50)) matches)]

            (is (> (count close-matches) 0)
                "Should find some peaks close to expected positions")

            (is (> (count close-matches) (* 0.5 (count expected-peaks)))
                (str "Should find majority of expected peaks: "
                     (count close-matches) " out of " (count expected-peaks)))))))))

(deftest test-noise-estimation-scaling
  (testing "Noise estimation scaling with different noise levels"
    (when (ensure-maldiquant-setup)
      (let [noise-levels [5.0 10.0 20.0 40.0]
            test-scenarios (map #(create-realistic-maldi-spectrum 300 :noise-level %) noise-levels)]

        (doseq [[noise-level test-data] (map vector noise-levels test-scenarios)]
          (let [test-spectrum (:spectrum test-data)
                mz-values (:mz-values test-spectrum)
                intensities (:intensities test-spectrum)

                ;; Clojure noise estimation
                clj-mad (peaks/estimate-noise intensities :mad)

                ;; R MALDIquant noise estimation
                _ (r "library(MALDIquant)")
                _ (r (str "mz <- " (to-r-vector mz-values)))
                _ (r (str "intensities <- " (to-r-vector intensities)))
                _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
                _ (r "r_noise <- estimateNoise(spectrum, method='MAD')")
                r-noise-data (-> (r "as.vector(r_noise)[-(1:length(mz))]") r/r->clj vec)
                r-mad (if (seq r-noise-data) (first r-noise-data) 0)]

            (is (> clj-mad 0)
                (str "Clojure MAD should be positive for noise level " noise-level))

            (is (> r-mad 0)
                (str "R MAD should be positive for noise level " noise-level))

            ;; Estimates should scale with noise level
            (when (> noise-level 10.0)
              (is (> clj-mad (* 0.3 noise-level))
                  (str "MAD should scale with noise level " noise-level ": " clj-mad)))

            ;; Cross-implementation consistency
            (let [ratio (/ clj-mad r-mad)]
              (is (and (> ratio 0.2) (< ratio 5.0))
                  (str "Implementations should be reasonably consistent for level "
                       noise-level ": " clj-mad " vs " r-mad)))))))))

(deftest test-statistical-reproducibility
  (testing "Statistical reproducibility across multiple runs"
    (when (ensure-maldiquant-setup)
      (let [n-replicates 3 ; Small number for fast testing
            test-spectra (repeatedly n-replicates #(create-realistic-maldi-spectrum 250))
            methods [:snip :tophat]]

        (doseq [method methods]
          (testing (str "method: " method)
            (let [correlations (map (fn [test-data]
                                      (let [test-spectrum (:spectrum test-data)
                                            mz-values (:mz-values test-spectrum)
                                            intensities (:intensities test-spectrum)

                                           ;; Clojure result
                                            clj-baseline (baseline/estimate-baseline test-spectrum method
                                                                                     :iterations 40 :half-window-size 20)

                                           ;; R result
                                            _ (r "library(MALDIquant)")
                                            _ (r (str "mz <- " (to-r-vector mz-values)))
                                            _ (r (str "intensities <- " (to-r-vector intensities)))
                                            _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
                                            _ (case method
                                                :snip (r "r_baseline <- estimateBaseline(spectrum, method='SNIP', iterations=40)")
                                                :tophat (r "r_baseline <- estimateBaseline(spectrum, method='TopHat', halfWindowSize=20)"))
                                            r-baseline (-> (r "as.vector(r_baseline)[-(1:length(mz))]") r/r->clj vec)

                                            clj-vec (if (dtype/reader? clj-baseline) (dtype/->vector clj-baseline) clj-baseline)]

                                        (calculate-correlation clj-vec r-baseline)))
                                    test-spectra)

                  mean-correlation (/ (reduce + correlations) (count correlations))]

              (is (> mean-correlation 0.4)
                  (str method " should have reasonable mean correlation: " mean-correlation))

              (is (every? #(> % 0.2) correlations)
                  (str method " should have consistent performance across replicates")))))))))

(deftest test-complete-workflow-comparison
  (testing "Complete MALDIquant workflow comparison"
    (when (ensure-maldiquant-setup)
      (let [test-data (create-realistic-maldi-spectrum 300
                                                       :peak-positions [1200 1600 2000 2400]
                                                       :peak-intensities [800 600 1000 500])
            test-spectrum (:spectrum test-data)
            mz-values (:mz-values test-spectrum)
            intensities (:intensities test-spectrum)]

        (testing "baseline + peak detection workflow"
          ;; Clojure workflow
          (let [clj-corrected (baseline/remove-baseline test-spectrum :snip :iterations 50)
                clj-peaks (peaks/detect-peaks clj-corrected
                                              :snr-threshold 3.0
                                              :window-half-size 8
                                              :noise-method :mad)
                clj-peak-count (peaks/peak-count clj-peaks)
                clj-peak-mz (peaks/get-peak-mz-values clj-peaks)

                ;; R workflow
                _ (r "library(MALDIquant)")
                _ (r (str "mz <- " (to-r-vector mz-values)))
                _ (r (str "intensities <- " (to-r-vector intensities)))
                _ (r "spectrum <- createMassSpectrum(mass=mz, intensity=intensities)")
                _ (r "corrected <- removeBaseline(spectrum, method='SNIP', iterations=50)")
                _ (r "r_peaks <- detectPeaks(corrected, method='MAD', SNR=3.0, halfWindowSize=8)")
                r-peak-mz (-> (r "mass(r_peaks)") r/r->clj vec)
                r-peak-count (count r-peak-mz)]

            (is (> clj-peak-count 0)
                "Clojure workflow should detect some peaks")

            (is (> r-peak-count 0)
                "R workflow should detect some peaks")

            ;; Peak counts should be reasonably similar
            (let [count-ratio (if (> r-peak-count 0) (/ (double clj-peak-count) (double r-peak-count)) 0)]
              (is (and (> count-ratio 0.3) (< count-ratio 3.0))
                  (str "Peak counts should be in reasonable range: " clj-peak-count " vs " r-peak-count)))

            ;; Should find peaks in expected regions
            (let [expected-peaks (:peak-positions test-data)
                  close-clj-matches (count (filter (fn [expected]
                                                     (some #(< (Math/abs (- % expected)) 50) clj-peak-mz))
                                                   expected-peaks))]
              (is (> close-clj-matches (* 0.5 (count expected-peaks)))
                  (str "Should find most expected peaks: " close-clj-matches " out of " (count expected-peaks))))))))))

;; ==============================================================================
;; Test Runner Functions  
;; ==============================================================================

(defn run-detailed-baseline-tests []
  "Run detailed baseline comparison tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running Detailed Baseline Tests ===")
    (clojure.test/run-test-var #'test-snip-parameter-sensitivity-detailed)
    (clojure.test/run-test-var #'test-baseline-methods-comprehensive-detailed)
    (println "=== Detailed baseline tests complete ===")))

(defn run-detailed-peak-tests []
  "Run detailed peak detection tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running Detailed Peak Detection Tests ===")
    (clojure.test/run-test-var #'test-peak-detection-parameter-analysis)
    (clojure.test/run-test-var #'test-noise-estimation-scaling)
    (println "=== Detailed peak tests complete ===")))

(defn run-workflow-tests []
  "Run complete workflow comparison tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running Workflow Comparison Tests ===")
    (clojure.test/run-test-var #'test-complete-workflow-comparison)
    (println "=== Workflow tests complete ===")))

(defn run-statistical-validation-tests []
  "Run statistical validation tests"
  (when (ensure-maldiquant-setup)
    (println "=== Running Statistical Validation Tests ===")
    (clojure.test/run-test-var #'test-statistical-reproducibility)
    (println "=== Statistical validation tests complete ===")))

(defn run-all-detailed-tests []
  "Run complete detailed MALDIquant comparison test suite"
  (when (ensure-maldiquant-setup)
    (println "=== Running Complete Detailed MALDIquant Test Suite ===")
    (run-detailed-baseline-tests)
    (run-detailed-peak-tests)
    (run-workflow-tests)
    (run-statistical-validation-tests)
    (println "=== Complete detailed test suite finished ==="))
  (when (not (ensure-maldiquant-setup))
    (println "MALDIquant R environment not available - skipping detailed tests")))
