(ns maldi-clj.alignment-test
  "Tests for peak alignment and binning algorithms"
  (:require [clojure.test :refer [deftest is testing]]
            [maldi-clj.alignment :as alignment]
            [maldi-clj.spectrum :as spectrum]
            [maldi-clj.peaks :as peaks]
            [tech.v3.datatype :as dtype]))

;; ==============================================================================
;; TEST UTILITIES
;; ==============================================================================

(defn create-test-mass-peaks
  "Create synthetic MassPeaks for testing alignment"
  []
  (let [;; Spectrum 1: peaks at m/z 100, 150, 200
        peaks1 {:id "test-spec-1"
                :peaks [{:mz 100.0 :intensity 1000 :snr 50.0 :noise 20.0 :index 50}
                        {:mz 150.0 :intensity 800 :snr 40.0 :noise 20.0 :index 500}
                        {:mz 200.0 :intensity 600 :snr 30.0 :noise 20.0 :index 1000}
                        {:mz 120.0 :intensity 300 :snr 15.0 :noise 20.0 :index 200}]
                :metadata {:test true}}

        ;; Spectrum 2: peaks at m/z 100.1, 150.2, 175 (slightly shifted + new peak)
        peaks2 {:id "test-spec-2"
                :peaks [{:mz 100.1 :intensity 900 :snr 45.0 :noise 20.0 :index 51}
                        {:mz 150.2 :intensity 700 :snr 35.0 :noise 20.0 :index 501}
                        {:mz 175.0 :intensity 500 :snr 25.0 :noise 20.0 :index 750}
                        {:mz 120.1 :intensity 280 :snr 14.0 :noise 20.0 :index 201}]
                :metadata {:test true}}

        ;; Spectrum 3: peaks at m/z 100.05, 175.1, 250 (partial overlap + new peak)
        peaks3 {:id "test-spec-3"
                :peaks [{:mz 100.05 :intensity 950 :snr 47.5 :noise 20.0 :index 50}
                        {:mz 175.1 :intensity 450 :snr 22.5 :noise 20.0 :index 751}
                        {:mz 250.0 :intensity 400 :snr 20.0 :noise 20.0 :index 1250}
                        {:mz 120.05 :intensity 290 :snr 14.5 :noise 20.0 :index 200}]
                :metadata {:test true}}]

    [peaks1 peaks2 peaks3]))

(defn create-test-spectra-with-peaks
  "Create test spectra and their detected peaks"
  []
  (let [;; Create synthetic spectra
        spec1 (spectrum/create-spectrum "spec-1"
                                        (range 90 260 0.1)
                                        (mapv (fn [mz]
                                                (+ 20 (* 5 (rand))
                                                   (cond
                                                     (< (Math/abs (- mz 100)) 0.5) 500
                                                     (< (Math/abs (- mz 150)) 0.5) 400
                                                     (< (Math/abs (- mz 200)) 0.5) 300
                                                     :else 0)))
                                              (range 90 260 0.1))
                                        {})

        spec2 (spectrum/create-spectrum "spec-2"
                                        (range 90 260 0.1)
                                        (mapv (fn [mz]
                                                (+ 20 (* 5 (rand))
                                                   (cond
                                                     (< (Math/abs (- mz 100.1)) 0.5) 450
                                                     (< (Math/abs (- mz 150.2)) 0.5) 350
                                                     (< (Math/abs (- mz 175)) 0.5) 250
                                                     :else 0)))
                                              (range 90 260 0.1))
                                        {})]

    [spec1 spec2]))

;; ==============================================================================
;; M/Z TOLERANCE TESTS
;; ==============================================================================

(deftest test-mz-tolerance-calculation
  (testing "M/Z tolerance calculation"
    (testing "Absolute tolerance"
      (let [[min-mz max-mz] (alignment/calculate-mz-tolerance 100.0 0.5 :absolute)]
        (is (= 99.5 min-mz) "Min m/z should be correct")
        (is (= 100.5 max-mz) "Max m/z should be correct")))

    (testing "PPM tolerance"
      (let [[min-mz max-mz] (alignment/calculate-mz-tolerance 1000.0 10 :ppm)]
        (is (= 999.99 min-mz) "Min m/z should be correct for PPM")
        (is (= 1000.01 max-mz) "Max m/z should be correct for PPM")))

    (testing "Invalid tolerance mode"
      (is (thrown? Exception (alignment/calculate-mz-tolerance 100.0 0.5 :invalid))
          "Should throw exception for invalid mode"))))

(deftest test-mz-within-tolerance
  (testing "M/Z tolerance checking"
    (testing "Within absolute tolerance"
      (is (alignment/mz-within-tolerance? 100.0 100.3 0.5 :absolute) "Should be within tolerance")
      (is (not (alignment/mz-within-tolerance? 100.0 101.0 0.5 :absolute)) "Should be outside tolerance"))

    (testing "Within PPM tolerance"
      (is (alignment/mz-within-tolerance? 1000.0 1000.005 10 :ppm) "Should be within PPM tolerance")
      (is (not (alignment/mz-within-tolerance? 1000.0 1000.05 10 :ppm)) "Should be outside PPM tolerance"))))

(deftest test-find-matching-peaks
  (testing "Finding matching peaks"
    (let [peak-list [{:mz 100.0} {:mz 150.0} {:mz 200.0}]
          matches (alignment/find-matching-peaks peak-list 100.2 0.5 :absolute)]
      (is (= [0] matches) "Should find first peak"))

    (let [mz-list [100.0 150.0 200.0]
          matches (alignment/find-matching-peaks mz-list 150.1 0.2 :absolute)]
      (is (= [1] matches) "Should find second peak in mz list"))))

;; ==============================================================================
;; PEAK BINNING TESTS
;; ==============================================================================

(deftest test-bin-peaks-strict
  (testing "Strict peak binning"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          binned (alignment/bin-peaks [peaks1 peaks2 peaks3]
                                      :tolerance 0.2
                                      :tolerance-mode :absolute
                                      :method :strict
                                      :min-frequency 2)]

      (testing "Basic structure"
        (is (map? binned) "Should return a map")
        (is (vector? (:bins binned)) "Should have bins vector")
        (is (= 3 (:total-spectra binned)) "Should track total spectra")
        (is (map? (:binning-parameters binned)) "Should have parameters"))

      (testing "Bins should contain expected peaks"
        (let [bins (:bins binned)
              bin-mz-centers (map :mz-center bins)]
          (is (> (count bins) 0) "Should find some bins")
          ;; Should find bins around m/z 100 and 120 (appear in multiple spectra)
          (is (some #(< (Math/abs (- % 100)) 1.0) bin-mz-centers) "Should find bin around m/z 100")
          (is (some #(< (Math/abs (- % 120)) 1.0) bin-mz-centers) "Should find bin around m/z 120")))

      (testing "Frequency filtering"
        (let [bins (:bins binned)]
          (is (every? #(>= (:frequency %) 2) bins) "All bins should meet min frequency"))))))

(deftest test-bin-peaks-reference
  (testing "Reference-based peak binning"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          binned (alignment/bin-peaks [peaks1 peaks2 peaks3]
                                      :tolerance 0.2
                                      :tolerance-mode :absolute
                                      :method :reference
                                      :reference-spectrum "test-spec-1"
                                      :min-frequency 1)]

      (testing "Reference mode structure"
        (is (= :reference (get-in binned [:binning-parameters :mode])) "Should use reference mode")
        (is (= "test-spec-1" (get-in binned [:binning-parameters :reference-spectrum])) "Should track reference"))

      (testing "Should align to reference peaks"
        (let [bins (:bins binned)
              reference-mz-values [100.0 150.0 200.0 120.0]]
          ;; Should find bins centered on reference m/z values
          (doseq [ref-mz reference-mz-values]
            (is (some #(= (:mz-center %) ref-mz) bins)
                (str "Should find bin centered on reference m/z " ref-mz))))))))

(deftest test-bin-peaks-parameters
  (testing "Peak binning parameter effects"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)]

      (testing "Tolerance effects"
        (let [strict-tolerance (alignment/bin-peaks [peaks1 peaks2 peaks3]
                                                    :tolerance 0.05
                                                    :method :strict
                                                    :min-frequency 2)
              loose-tolerance (alignment/bin-peaks [peaks1 peaks2 peaks3]
                                                   :tolerance 0.5
                                                   :method :strict
                                                   :min-frequency 2)]
          (is (<= (count (:bins strict-tolerance)) (count (:bins loose-tolerance)))
              "Stricter tolerance should produce same or fewer bins")))

      (testing "Minimum frequency effects"
        (let [low-freq (alignment/bin-peaks [peaks1 peaks2 peaks3]
                                            :min-frequency 1)
              high-freq (alignment/bin-peaks [peaks1 peaks2 peaks3]
                                             :min-frequency 3)]
          (is (>= (count (:bins low-freq)) (count (:bins high-freq)))
              "Lower frequency threshold should produce more bins"))))))

;; ==============================================================================
;; PEAK FILTERING TESTS
;; ==============================================================================

(deftest test-filter-peaks-by-frequency
  (testing "Peak filtering by frequency"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          binned (alignment/bin-peaks [peaks1 peaks2 peaks3] :min-frequency 1)
          filtered (alignment/filter-peaks-by-frequency binned 3 :frequency-mode :absolute)]

      (testing "Frequency filtering"
        (is (<= (count (:bins filtered)) (count (:bins binned)))
            "Filtering should reduce or maintain bin count")
        (is (every? #(>= (:frequency %) 3) (:bins filtered))
            "All remaining bins should meet frequency threshold"))

      (testing "Relative frequency mode"
        (let [relative-filtered (alignment/filter-peaks-by-frequency binned 0.8 :frequency-mode :relative)]
          (is (every? #(>= (:frequency %) 2) (:bins relative-filtered))
              "Should convert relative frequency correctly"))))))

(deftest test-filter-peaks-by-intensity
  (testing "Peak filtering by intensity"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          binned (alignment/bin-peaks [peaks1 peaks2 peaks3] :min-frequency 1)
          filtered (alignment/filter-peaks-by-intensity binned 500 :intensity-mode :any)]

      (testing "Intensity filtering with :any mode"
        (is (<= (count (:bins filtered)) (count (:bins binned)))
            "Filtering should reduce or maintain bin count")
        (is (every? (fn [bin]
                      (some #(>= (:intensity %) 500) (:peaks bin)))
                    (:bins filtered))
            "Each bin should have at least one peak above threshold"))

      (testing "Intensity filtering with :all mode"
        (let [all-filtered (alignment/filter-peaks-by-intensity binned 100 :intensity-mode :all)]
          (is (every? (fn [bin]
                        (every? #(>= (:intensity %) 100) (:peaks bin)))
                      (:bins all-filtered))
              "All peaks in each bin should be above threshold"))))))

;; ==============================================================================
;; UTILITY FUNCTION TESTS
;; ==============================================================================

(deftest test-bin-statistics
  (testing "Bin statistics calculation"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          binned (alignment/bin-peaks [peaks1 peaks2 peaks3] :min-frequency 1)
          stats (alignment/bin-statistics binned)]

      (testing "Statistics structure"
        (is (map? stats) "Should return statistics map")
        (is (number? (:total-bins stats)) "Should have total bins count")
        (is (number? (:total-peaks stats)) "Should have total peaks count")
        (is (map? (:frequency-stats stats)) "Should have frequency statistics")
        (is (map? (:mz-range stats)) "Should have m/z range statistics"))

      (testing "Statistics values"
        (is (= (:total-bins stats) (count (:bins binned))) "Total bins should match")
        (is (> (:total-peaks stats) 0) "Should have positive peak count")))))

(deftest test-merge-mass-peaks
  (testing "Mass peaks merging"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          merged (alignment/merge-mass-peaks [peaks1 peaks2] :merge-strategy :union)]

      (testing "Union merge"
        (is (string? (:id merged)) "Should have merged ID")
        (is (vector? (:peaks merged)) "Should have peaks vector")
        (is (= (count (:peaks merged))
               (+ (count (:peaks peaks1)) (count (:peaks peaks2))))
            "Should combine all peaks")
        (is (= :union (get-in merged [:metadata :merge-strategy])) "Should track merge strategy")))))

(deftest test-create-intensity-matrix
  (testing "Intensity matrix creation"
    (let [[peaks1 peaks2 peaks3] (create-test-mass-peaks)
          binned (alignment/bin-peaks [peaks1 peaks2 peaks3] :min-frequency 1)
          matrix (alignment/create-intensity-matrix binned :missing-value-strategy :zero)]

      (testing "Matrix structure"
        (is (map? matrix) "Should return matrix map")
        (is (vector? (:spectrum-ids matrix)) "Should have spectrum IDs")
        (is (vector? (:column-names matrix)) "Should have column names")
        (is (vector? (:matrix-data matrix)) "Should have matrix data")
        (is (vector? (:dimensions matrix)) "Should have dimensions"))

      (testing "Matrix dimensions"
        (let [[rows cols] (:dimensions matrix)]
          (is (= rows 3) "Should have 3 rows (spectra)")
          (is (= cols (count (:bins binned))) "Should have correct number of columns")
          (is (= (count (:matrix-data matrix)) rows) "Matrix data should match row count")))

      (testing "Missing value handling"
        (let [zero-matrix (alignment/create-intensity-matrix binned :missing-value-strategy :zero)
              na-matrix (alignment/create-intensity-matrix binned :missing-value-strategy :na)]
          (is (every? (fn [row] (every? #(or (number? %) (= % ##NaN)) (:intensities row)))
                      (:matrix-data zero-matrix))
              "Zero strategy should produce numbers")
          (is (some (fn [row] (some #(Double/isNaN %) (:intensities row)))
                    (:matrix-data na-matrix))
              "NA strategy should produce some NaN values"))))))

;; ==============================================================================
;; INTEGRATION TESTS
;; ==============================================================================

(deftest test-full-alignment-pipeline
  (testing "Complete alignment pipeline"
    (let [[spec1 spec2] (create-test-spectra-with-peaks)
          ;; Detect peaks
          peaks1 (peaks/detect-peaks spec1 :snr-threshold 5.0 :baseline-corrected? true)
          peaks2 (peaks/detect-peaks spec2 :snr-threshold 5.0 :baseline-corrected? true)
          ;; Bin peaks
          binned (alignment/bin-peaks [peaks1 peaks2]
                                      :tolerance 0.5
                                      :method :strict
                                      :min-frequency 1)
          ;; Filter by frequency
          filtered (alignment/filter-peaks-by-frequency binned 2)
          ;; Create intensity matrix
          matrix (alignment/create-intensity-matrix filtered)]

      (testing "Pipeline produces valid results"
        (is (> (peaks/peak-count peaks1) 0) "Should detect peaks in first spectrum")
        (is (> (peaks/peak-count peaks2) 0) "Should detect peaks in second spectrum")
        (is (> (count (:bins binned)) 0) "Should create bins")
        (is (<= (count (:bins filtered)) (count (:bins binned))) "Filtering should work")
        (is (map? matrix) "Should create intensity matrix"))

      (testing "Expected peaks should be aligned"
        (let [bins (:bins filtered)]
          ;; Should find aligned peaks around the major synthetic peaks
          (is (some (fn [bin]
                      (and (< (Math/abs (- (:mz-center bin) 100)) 5)
                           (>= (:frequency bin) 2)))
                    bins)
              "Should find aligned peaks around m/z 100"))))))

;; ==============================================================================
;; ERROR HANDLING TESTS
;; ==============================================================================

(deftest test-error-handling
  (testing "Error handling"
    (testing "Empty mass peaks collection"
      (is (thrown? Exception (alignment/bin-peaks []))
          "Should throw exception for empty collection"))

    (testing "Invalid binning method"
      (let [[peaks1 peaks2] (create-test-mass-peaks)]
        (is (thrown? Exception (alignment/bin-peaks [peaks1 peaks2] :method :invalid))
            "Should throw exception for invalid method")))

    (testing "Reference spectrum not found"
      (let [[peaks1 peaks2] (create-test-mass-peaks)]
        (is (thrown? Exception (alignment/bin-peaks [peaks1 peaks2]
                                                    :method :reference
                                                    :reference-spectrum "nonexistent"))
            "Should throw exception for missing reference")))))

;; ==============================================================================
;; TEST RUNNER
;; ==============================================================================

(defn run-alignment-tests
  "Run all alignment tests"
  []
  (clojure.test/run-tests 'maldi-clj.alignment-test))
