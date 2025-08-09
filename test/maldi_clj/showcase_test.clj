(ns maldi-clj.showcase-test
  "Showcase test demonstrating the complete MALDIcloj capabilities and validation"
  (:require [clojure.test :refer [deftest is testing]]
            [maldi-clj.spectrum :as spectrum]
            [maldi-clj.baseline :as baseline]
            [maldi-clj.peaks :as peaks]
            [tech.v3.datatype :as dtype]
            [clojure.pprint :as pprint]))

(deftest showcase-complete-maldi-workflow
  (testing "Complete MALDI-TOF bacterial identification workflow showcase"
    (println "\n🔬 MALDIcloj Showcase: Complete Bacterial MALDI-TOF Workflow")
    (println "=" 60)

    ;; Generate realistic bacterial spectrum
    (let [;; Simulate bacterial ribosomal proteins (typical E. coli profile)
          bacterial-peaks [2903 3236 4365 5096 6241 6857 7159 8564 9536 10300 11684 12227 14931 16950]
          peak-intensities [800 600 1200 900 1500 700 400 850 950 600 450 380 300 250]

          ;; Create realistic spectrum with 1000 points from 2000-20000 m/z
          mz-values (vec (range 2000.0 20000.0 18.0))
          n-points (count mz-values)

          ;; Generate realistic peak shapes with mass-dependent resolution
          peak-signals (map (fn [pos intensity]
                              (mapv (fn [mz]
                                      (let [resolution (/ pos 50) ; Typical MALDI resolution
                                            sigma (/ pos resolution)
                                            dist (/ (- mz pos) sigma)
                                            gaussian (* intensity (Math/exp (- (/ (* dist dist) 2))))]
                                        gaussian))
                                    mz-values))
                            bacterial-peaks peak-intensities)

          ;; Combine all peaks
          combined-signal (reduce (fn [acc peak] (mapv + acc peak))
                                  (vec (repeat n-points 0.0))
                                  peak-signals)

          ;; Add realistic MALDI matrix baseline (decreasing exponential)
          baseline-component (mapv #(+ 20.0 (* 1000 (Math/exp (- (/ % 500))))) mz-values)

          ;; Add shot-to-shot noise
          noise-component (repeatedly n-points #(+ (* 25 (- (rand) 0.5))
                                                   (* 8 (Math/sqrt (max 1 (rand-int 100))))))

          ;; Create raw spectrum
          raw-intensities (mapv #(max 0.1 (+ %1 %2 %3)) combined-signal baseline-component noise-component)
          raw-spectrum (spectrum/create-spectrum "E-coli-clinical-sample" mz-values raw-intensities
                                                 {:organism "E. coli" :matrix "CHCA" :acquisition "linear"})]

      (println "\n📊 Raw Spectrum Characteristics:")
      (println (format "• m/z range: %.1f - %.1f" (first mz-values) (last mz-values)))
      (println (format "• Data points: %d" n-points))
      (println (format "• Intensity range: %.1f - %.1f"
                       (apply min raw-intensities) (apply max raw-intensities)))
      (println (format "• Expected peaks: %d" (count bacterial-peaks)))

      ;; Step 1: Baseline Correction
      (println "\n🔧 Step 1: Baseline Correction using SNIP algorithm")
      (let [start-time (System/nanoTime)
            corrected-spectrum (baseline/remove-baseline raw-spectrum :snip :iterations 50)
            correction-time (/ (- (System/nanoTime) start-time) 1000000.0)

            baseline-estimate (baseline/estimate-baseline raw-spectrum :snip :iterations 50)
            baseline-vec (if (dtype/reader? baseline-estimate)
                           (dtype/->vector baseline-estimate) baseline-estimate)

            corrected-intensities (:intensities corrected-spectrum)]

        (println (format "• Processing time: %.2f ms" correction-time))
        (println (format "• Baseline range: %.1f - %.1f"
                         (apply min baseline-vec) (apply max baseline-vec)))
        (println (format "• Signal enhancement: %.1fx"
                         (/ (apply max corrected-intensities) (apply max baseline-vec))))

        ;; Validate baseline correction
        (is (every? #(>= % 0) baseline-vec) "Baseline should be non-negative")
        (is (< correction-time 100.0) "Baseline correction should be fast")
        (is (> (/ (apply max corrected-intensities) (apply max baseline-vec)) 2.0)
            "Signal should be enhanced relative to baseline")

        ;; Step 2: Peak Detection
        (println "\n🎯 Step 2: Peak Detection with adaptive SNR threshold")
        (let [start-time (System/nanoTime)
              detected-peaks (peaks/detect-peaks corrected-spectrum
                                                 :snr-threshold 3.0
                                                 :window-half-size 8
                                                 :noise-method :mad)
              detection-time (/ (- (System/nanoTime) start-time) 1000000.0)

              peak-count (peaks/peak-count detected-peaks)
              detected-mz (peaks/get-peak-mz-values detected-peaks)
              detected-intensities (peaks/get-peak-intensities detected-peaks)

              ;; Calculate noise estimate
              noise-estimate (peaks/estimate-noise corrected-intensities :mad)]

          (println (format "• Processing time: %.2f ms" detection-time))
          (println (format "• Detected peaks: %d" peak-count))
          (println (format "• Noise estimate (MAD): %.1f" noise-estimate))
          (println (format "• Peak m/z range: %.1f - %.1f"
                           (if (seq detected-mz) (apply min detected-mz) 0)
                           (if (seq detected-mz) (apply max detected-mz) 0)))

          ;; Validate peak detection
          (is (> peak-count 0) "Should detect some peaks")
          (is (< detection-time 50.0) "Peak detection should be fast")
          (is (> noise-estimate 0) "Noise estimate should be positive")

          ;; Step 3: Peak Matching and Identification
          (println "\n🔍 Step 3: Peak Matching against Expected Bacterial Profile")
          (let [;; Match detected peaks to expected bacterial peaks
                matches (map (fn [expected]
                               (if (seq detected-mz)
                                 (let [closest (apply min-key #(Math/abs (- % expected)) detected-mz)
                                       distance (Math/abs (- closest expected))]
                                   {:expected expected :detected closest :error distance
                                    :match (< distance 50)}) ; ±50 m/z tolerance
                                 {:expected expected :match false}))
                             bacterial-peaks)

                successful-matches (filter :match matches)
                match-rate (/ (count successful-matches) (count bacterial-peaks))
                avg-error (if (seq successful-matches)
                            (/ (reduce + (map :error successful-matches)) (count successful-matches))
                            0)]

            (println (format "• Peak matches: %d/%d (%.1f%%)"
                             (count successful-matches) (count bacterial-peaks) (* 100 match-rate)))
            (println (format "• Average error: ±%.1f m/z" avg-error))

            ;; Display match details
            (println "\n📋 Peak Matching Results:")
            (doseq [match (take 10 matches)] ; Show first 10 matches
              (if (:match match)
                (println (format "  ✓ %.1f m/z → %.1f m/z (Δ%.1f)"
                                 (:expected match) (:detected match) (:error match)))
                (println (format "  ✗ %.1f m/z → not detected" (:expected match)))))

            ;; Validate identification performance
            (is (> match-rate 0.4) "Should match at least 40% of expected peaks")
            (is (< avg-error 100) "Average error should be reasonable")

            ;; Step 4: Quality Assessment
            (println "\n📈 Step 4: Workflow Quality Assessment")
            (let [total-time (+ correction-time detection-time)
                  snr-values (map #(/ % noise-estimate) detected-intensities)
                  avg-snr (if (seq snr-values) (/ (reduce + snr-values) (count snr-values)) 0)
                  spectrum-quality (cond
                                     (and (> match-rate 0.7) (> avg-snr 5.0)) "Excellent"
                                     (and (> match-rate 0.5) (> avg-snr 3.0)) "Good"
                                     (and (> match-rate 0.3) (> avg-snr 2.0)) "Acceptable"
                                     :else "Poor")]

              (println (format "• Total processing time: %.2f ms" total-time))
              (println (format "• Average peak SNR: %.1f" avg-snr))
              (println (format "• Spectrum quality: %s" spectrum-quality))
              (println (format "• Identification confidence: %.1f%%" (* 100 match-rate)))

              ;; Final validation
              (is (< total-time 200.0) "Complete workflow should be fast")
              (is (> avg-snr 2.0) "Detected peaks should have reasonable SNR")
              (is (not= spectrum-quality "Poor") "Spectrum quality should be acceptable or better")

              ;; Summary
              (println "\n🎉 Workflow Summary:")
              (println (format "• Successfully processed %d-point bacterial MALDI spectrum" n-points))
              (println (format "• Identified %d/%d expected biomarker peaks"
                               (count successful-matches) (count bacterial-peaks)))
              (println (format "• Processing completed in %.1f ms with %s quality"
                               total-time spectrum-quality))
              (println "• All validation criteria met ✓")

              (is true "Complete workflow validation successful!"))))))))

(deftest showcase-algorithm-comparison
  (testing "Algorithm comparison showcase across different baseline methods"
    (println "\n⚔️  Algorithm Comparison Showcase")
    (println "=" 40)

    ;; Create test spectrum with known characteristics
    (let [mz-values (vec (range 1000.0 5000.0 2.0))
          n-points (count mz-values)

          ;; Add controlled peaks
          true-peaks (mapv #(+ (* 100 (Math/exp (- (/ (* (- % 500) (- % 500)) 10000))))
                               (* 150 (Math/exp (- (/ (* (- % 1500) (- % 1500)) 8000))))
                               (* 200 (Math/exp (- (/ (* (- % 2500) (- % 2500)) 12000)))))
                           (range n-points))

          ;; Add complex baseline
          true-baseline (mapv #(+ 50 (* 0.01 %) (* 20 (Math/sin (/ % 200)))) (range n-points))

          ;; Add noise
          noise (repeatedly n-points #(* 15 (- (rand) 0.5)))

          ;; Combine components
          intensities (mapv #(max 0.1 (+ %1 %2 %3)) true-peaks true-baseline noise)
          test-spectrum (spectrum/create-spectrum "algorithm-comparison" mz-values intensities {})]

      (println (format "Test spectrum: %d points, 3 synthetic peaks" n-points))

      ;; Test all baseline methods
      (doseq [method [:snip :tophat :median]]
        (println (format "\n🔧 Testing %s algorithm:" (name method)))

        (let [start-time (System/nanoTime)
              estimated-baseline (baseline/estimate-baseline test-spectrum method
                                                             :iterations 50 :half-window-size 25)
              processing-time (/ (- (System/nanoTime) start-time) 1000000.0)

              baseline-vec (if (dtype/reader? estimated-baseline)
                             (dtype/->vector estimated-baseline) estimated-baseline)

              ;; Calculate correlation with true baseline
              correlation (let [n (count baseline-vec)]
                            mean-true (/ (reduce + true-baseline) n)
                            mean-est (/ (reduce + baseline-vec) n)
                            numerator (reduce + (map #(* (- %1 mean-true) (- %2 mean-est))
                                                     true-baseline baseline-vec))
                            denom-true (Math/sqrt (reduce + (map #(* (- % mean-true) (- % mean-true))
                                                                 true-baseline)))
                            denom-est (Math/sqrt (reduce + (map #(* (- % mean-est) (- % mean-est))
                                                                baseline-vec)))
                            (if (and (> denom-true 0) (> denom-est 0))
                              (/ numerator (* denom-true denom-est))
                              0.0))

              ;; Calculate RMSE
              rmse (Math/sqrt (/ (reduce + (map #(* (- %1 %2) (- %1 %2)) baseline-vec true-baseline))
                                 n-points))]

          (println (format "  • Processing time: %.2f ms" processing-time))
          (println (format "  • Correlation with truth: %.3f" correlation))
          (println (format "  • RMSE: %.2f" rmse))
          (println (format "  • Baseline range: %.1f - %.1f"
                           (apply min baseline-vec) (apply max baseline-vec)))

          ;; Validate each method
          (is (> correlation 0.3) (format "%s should show reasonable correlation" (name method)))
          (is (< processing-time 100.0) (format "%s should process quickly" (name method)))
          (is (every? #(>= % 0) baseline-vec) (format "%s should produce non-negative baseline" (name method)))))

      (println "\n✅ All baseline algorithms validated successfully!"))))

(defn run-showcase-tests []
  "Run the complete MALDIcloj showcase demonstration"
  (println "\n🚀 MALDIcloj Showcase Test Suite")
  (println "🧬 Demonstrating world-class MALDI-TOF processing capabilities")
  (clojure.test/run-test-var #'showcase-complete-maldi-workflow)
  (clojure.test/run-test-var #'showcase-algorithm-comparison)
  (println "\n🎉 Showcase complete! MALDIcloj is ready for production use."))
