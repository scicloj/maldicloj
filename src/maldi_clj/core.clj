(ns maldi-clj.core
  "Core API for MALDI-TOF mass spectrometry data analysis.
   
   This namespace provides the main entry points for MaldiCloj functionality,
   offering simplified high-level functions for common workflows."
  (:require [maldi-clj.spectrum :as spectrum]
            [maldi-clj.baseline :as baseline]
            [maldi-clj.peaks :as peaks]
            [maldi-clj.smoothing :as smoothing]
            [maldi-clj.preprocessing :as preprocessing]
            [maldi-clj.alignment :as alignment]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [tech.v3.dataset :as ds]))

;; =============================================================================
;; Core Spectrum Creation and Validation
;; =============================================================================

(defn create-spectrum
  "Create a MALDI-TOF spectrum with automatic performance optimization.
   
   Automatically chooses optimal data types and memory layouts:
   - Uses float32 for intensities (2x memory savings, minimal precision loss)
   - Uses native buffers for very large datasets (>20k points) for memory efficiency
   - Provides explicit control through options for advanced users
   
   Options:
   - :validate? - Validate the created spectrum (default: true)
   - :optimize? - Use optimized data structures (default: true, if false creates legacy spectrum)
   - :add-performance-metadata? - Add performance info to metadata (default: true)
   - :mz-dtype - Data type for m/z values (:float32 or :float64)
   - :intensity-dtype - Data type for intensities (:float32 or :float64)
   - :use-native-buffers? - Use native memory buffers (auto-configured if not specified)"
  [id mz-values intensities & [metadata & options]]
  (let [metadata (or metadata {})
        n-points (count mz-values)

        ;; Parse options
        {:keys [validate? optimize? add-performance-metadata? mz-dtype intensity-dtype use-native-buffers?]
         :or {validate? true optimize? true add-performance-metadata? true}} options]

    ;; Validation for invalid inputs
    (when validate?
      (when (empty? mz-values)
        (throw (ex-info "Cannot create spectrum with empty m/z values"
                        {:mz-count (count mz-values) :intensity-count (count intensities)})))
      (when (not= (count mz-values) (count intensities))
        (throw (ex-info "m/z values and intensities must have the same length"
                        {:mz-count (count mz-values) :intensity-count (count intensities)}))))

    ;; Handle legacy spectrum creation
    (if (not optimize?)
      (spectrum/create-spectrum-legacy id mz-values intensities metadata)

      ;; Optimized spectrum creation
      (let [;; Conservative auto-configuration - prioritize reliability over speed
            auto-config (cond
                          ;; Very large spectra - use native buffers for memory efficiency
                          (>= n-points 100000)
                          {:mz-dtype :float32 ; Trade precision for memory on huge datasets
                           :intensity-dtype :float32 ; Always use float32 for intensities  
                           :use-native-buffers? true} ; Native memory for memory efficiency

                          ;; Large spectra - balance memory and precision
                          (>= n-points 20000)
                          {:mz-dtype :float64 ; Keep precision for m/z values
                           :intensity-dtype :float32 ; Optimize intensities
                           :use-native-buffers? true} ; Native memory for better memory management

                          ;; Medium/small spectra - prioritize simplicity and speed
                          :else
                          {:mz-dtype :float64 ; Full precision
                           :intensity-dtype :float32 ; Still optimize intensities (main benefit)
                           :use-native-buffers? false}) ; JVM heap is faster for smaller datasets

            ;; Override auto-config with explicit options
            final-mz-dtype (or mz-dtype (:mz-dtype auto-config))
            final-intensity-dtype (or intensity-dtype (:intensity-dtype auto-config))
            final-use-native? (if (some? use-native-buffers?)
                                use-native-buffers?
                                (:use-native-buffers? auto-config))

            container-type (if final-use-native? :native-buffer :jvm-heap)

            spectrum (spectrum/create-spectrum id mz-values intensities metadata
                                               :mz-dtype final-mz-dtype
                                               :intensity-dtype final-intensity-dtype
                                               :container-type container-type)

            ;; Conditionally add performance metadata for introspection
            enhanced-spectrum (if add-performance-metadata?
                                (assoc-in spectrum [:metadata :performance]
                                          {:n-points n-points
                                           :mz-dtype final-mz-dtype
                                           :intensity-dtype final-intensity-dtype
                                           :native-buffers? final-use-native?
                                           :container-type container-type
                                           :auto-optimized? true
                                           :optimization-rationale
                                           (cond
                                             (>= n-points 100000) "Very large dataset - prioritize memory efficiency"
                                             (>= n-points 20000) "Large dataset - balance memory and precision"
                                             :else "Small/medium dataset - prioritize speed and simplicity")})
                                spectrum)]

        (if validate?
          (if (spectrum/valid-spectrum? enhanced-spectrum)
            enhanced-spectrum
            (throw (ex-info "Invalid spectrum created"
                            {:errors (spectrum/explain-spectrum-errors enhanced-spectrum)
                             :spectrum enhanced-spectrum})))
          enhanced-spectrum)))))

(defn valid-spectrum?
  "Check if spectrum is valid according to MaldiCloj standards."
  [spectrum]
  (spectrum/valid-spectrum? spectrum))

(defn spectrum-info
  "Get comprehensive information about a spectrum including performance characteristics."
  [spectrum]
  (let [;; Extract base metadata without performance info
        base-metadata (dissoc (:metadata spectrum) :performance)

        base-info {:id (:id spectrum)
                   :metadata base-metadata
                   :length (spectrum/spectrum-length spectrum)
                   :mz-range (spectrum/get-mz-range spectrum)
                   :intensity-stats (spectrum/get-intensity-stats spectrum)
                   :total-ion-current (preprocessing/total-ion-current spectrum)
                   :is-empty? (preprocessing/is-empty-spectrum? spectrum)
                   :dtype-optimized? (spectrum/dtype-spectrum? spectrum)}

        ;; Add performance information if available
        perf-info (get-in spectrum [:metadata :performance])

        memory-estimate (when perf-info
                          (let [n-points (:n-points perf-info)
                                mz-bytes (* n-points (if (= (:mz-dtype perf-info) :float32) 4 8))
                                intensity-bytes (* n-points 4)]
                            {:estimated-memory-mb (/ (+ mz-bytes intensity-bytes) 1024.0 1024.0)
                             :mz-memory-mb (/ mz-bytes 1024.0 1024.0)
                             :intensity-memory-mb (/ intensity-bytes 1024.0 1024.0)}))]

    (cond-> base-info
      perf-info (assoc :performance perf-info)
      memory-estimate (assoc :memory-usage memory-estimate))))

;; =============================================================================
;; High-Level Preprocessing Functions
;; =============================================================================

(defn clinical-preprocess
  "Apply standard clinical preprocessing pipeline (Weis et al. 2020)."
  [spectrum & {:keys [trim-range savgol-window savgol-order snip-iterations tic-normalize?]
               :or {trim-range [2000 20000] savgol-window 21 savgol-order 3
                    snip-iterations 20 tic-normalize? true}}]
  (when (preprocessing/is-empty-spectrum? spectrum)
    (throw (ex-info "Cannot preprocess empty spectrum" {:spectrum spectrum})))

  (let [mz-range (spectrum/get-mz-range spectrum)]
    (when (and trim-range
               (or (> (first trim-range) (second mz-range))
                   (< (second trim-range) (first mz-range))))
      (throw (ex-info "Trim range does not overlap with spectrum mass range"
                      {:trim-range trim-range
                       :spectrum-range mz-range
                       :spectrum-id (:id spectrum)}))))

  (-> spectrum
      ;; Step 1: Variance stabilization + Smoothing
      (smoothing/clinical-preprocessing {:savgol-window savgol-window
                                         :savgol-order savgol-order})
      ;; Step 2: Baseline correction
      (baseline/remove-baseline :snip :iterations snip-iterations)
      ;; Step 3: TIC normalization
      (cond-> tic-normalize?
        (preprocessing/preprocess-spectrum :tic-normalize? true))
      ;; Step 4: Mass range trimming
      (preprocessing/trim-spectrum (first trim-range) (second trim-range))))

(defn basic-preprocess
  "Apply basic preprocessing for general use."
  [spectrum & {:keys [smooth? smooth-method baseline-correct? baseline-method
                      normalize? normalize-method]
               :or {smooth? true smooth-method :savitzky-golay
                    baseline-correct? true baseline-method :snip
                    normalize? true normalize-method :tic}}]
  (cond-> spectrum
    smooth?
    ((fn [s]
       (case smooth-method
         :savitzky-golay (smoothing/savitzky-golay s)
         :moving-average (smoothing/moving-average s 11)
         s)))

    baseline-correct?
    (baseline/remove-baseline baseline-method)

    normalize?
    (preprocessing/preprocess-spectrum :tic-normalize? true)))

;; =============================================================================
;; Peak Detection and Analysis
;; =============================================================================

(defn detect-peaks
  "Detect peaks in a spectrum with sensible defaults."
  [spectrum & {:keys [snr-threshold window-size noise-method preprocess?]
               :or {snr-threshold 3.0 window-size 10 noise-method :mad preprocess? true}}]
  (let [processed-spectrum (if preprocess?
                             (try
                               (basic-preprocess spectrum)
                               (catch Exception e
                                 ;; If preprocessing fails, use the original spectrum
                                 (println "Warning: Preprocessing failed, using original spectrum for peak detection")
                                 spectrum))
                             spectrum)
        peak-results (peaks/detect-peaks processed-spectrum
                                         :snr-threshold snr-threshold
                                         :window-half-size window-size
                                         :noise-method noise-method)]
    ;; Return just the vector of peaks for simpler API
    (:peaks peak-results)))

(defn peak-statistics
  "Get comprehensive peak statistics from a peaks vector."
  [peaks]
  (when (seq peaks)
    (let [intensities (map :intensity peaks)
          snr-values (map :snr peaks)]
      {:count (count peaks)
       :intensity-stats {:min (apply min intensities)
                         :max (apply max intensities)
                         :mean (/ (reduce + intensities) (count intensities))}
       :snr-stats {:min (apply min snr-values)
                   :max (apply max snr-values)
                   :mean (/ (reduce + snr-values) (count snr-values))}})))

;; =============================================================================
;; Spectrum Comparison and Alignment
;; =============================================================================

(defn compare-spectra
  "Compare two spectra and return similarity metrics.
   
   Note: Full implementation pending alignment namespace completion."
  [spectrum1 spectrum2 & {:keys [preprocess? align? tolerance-ppm]
                          :or {preprocess? true align? true tolerance-ppm 500}}]
  (let [spec1 (if preprocess? (basic-preprocess spectrum1) spectrum1)
        spec2 (if preprocess? (basic-preprocess spectrum2) spectrum2)]
    {:status "placeholder"
     :message "Spectrum comparison not yet fully implemented"
     :spectrum1-info (spectrum-info spec1)
     :spectrum2-info (spectrum-info spec2)}))

;; =============================================================================
;; Batch Processing Functions  
;; =============================================================================

(defn process-batch
  "Process multiple spectra with the same parameters."
  [spectra processing-fn & {:keys [parallel? progress?]
                            :or {parallel? true progress? false}}]
  (let [process-fn (if parallel? pmap map)
        results (process-fn processing-fn spectra)]
    (if progress?
      (do (println (str "Processed " (count results) " spectra"))
          (vec results))
      (vec results))))

(defn batch-preprocess
  "Preprocess multiple spectra with optimal performance settings.
   
   Automatically configures batch size and parallelization based on:
   - Number of spectra
   - System memory and CPU cores
   - Spectrum sizes"
  [spectra & {:keys [method batch-size parallel?] :or {method :basic} :as options}]
  (let [n-spectra (count spectra)

        ;; Auto-configure batch processing based on system resources
        runtime (Runtime/getRuntime)
        max-memory-mb (/ (.maxMemory runtime) 1024 1024)
        cores (.availableProcessors runtime)

        ;; Smart defaults for batch processing
        auto-batch-size (cond
                          (>= max-memory-mb 16000) 100 ; High-memory systems
                          (>= max-memory-mb 8000) 50 ; Medium-memory systems  
                          :else 25) ; Conservative for low-memory

        final-batch-size (or batch-size auto-batch-size)
        final-parallel? (if (some? parallel?) parallel? (> n-spectra 10))

        process-fn (case method
                     :clinical #(apply clinical-preprocess % (mapcat identity options))
                     :basic #(apply basic-preprocess % (mapcat identity options)))

        ;; Use optimized batch processing
        processing-fn (if final-parallel? pmap map)]

    ;; Process in batches to control memory usage
    (->> spectra
         (partition-all final-batch-size)
         (processing-fn (fn [batch]
                          (mapv process-fn batch)))
         (apply concat)
         (vec))))

(defn batch-detect-peaks
  "Detect peaks in multiple spectra."
  [spectra & options]
  (let [detect-fn #(apply detect-peaks % (mapcat identity options))]
    (process-batch spectra detect-fn)))

;; =============================================================================
;; Data Export and Conversion
;; =============================================================================

(defn spectra->dataset
  "Convert spectra to a dataset for analysis.
   
   Two modes:
   - include-peaks? false: One row per spectrum with summary info
   - include-peaks? true: One row per spectrum with peak information included"
  [spectra & {:keys [include-peaks? peak-options]
              :or {include-peaks? false peak-options {}}}]
  (let [rows (map (fn [spectrum]
                    (let [info (spectrum-info spectrum)
                          base-row {:spectrum-id (:id spectrum)
                                    :length (:length info)
                                    :mz-range-start (first (:mz-range info))
                                    :mz-range-end (second (:mz-range info))
                                    :total-ion-current (:total-ion-current info)
                                    :is-empty (:is-empty? info)
                                    :dtype-optimized (:dtype-optimized? info)}
                          metadata-row (dissoc (:metadata spectrum) :performance)]
                      (if include-peaks?
                        (let [peaks (apply detect-peaks spectrum
                                           (mapcat identity peak-options))
                              peak-stats (peak-statistics peaks)]
                          (merge base-row
                                 metadata-row
                                 {:peak-count (count peaks)
                                  :peak-intensity-max (get-in peak-stats [:intensity-stats :max] 0)
                                  :peak-intensity-mean (get-in peak-stats [:intensity-stats :mean] 0)
                                  :peak-snr-max (get-in peak-stats [:snr-stats :max] 0)
                                  :peak-snr-mean (get-in peak-stats [:snr-stats :mean] 0)}))
                        (merge base-row metadata-row))))
                  spectra)]
    (ds/->>dataset rows)))

(defn export-results
  "Export analysis results to various formats."
  [data format filepath]
  (case format
    :csv (if (ds/dataset? data)
           (ds/write! data filepath)
           (throw (ex-info "CSV export requires dataset format" {:data data})))
    :edn (spit filepath (pr-str data))
    :json (throw (ex-info "JSON export not yet implemented" {:format format})))
  (str "Data exported to " filepath " in " (name format) " format"))

;; =============================================================================
;; Utility and Convenience Functions
;; =============================================================================

(defn benchmark-processing
  "Benchmark processing performance on test data."
  [n-spectra n-points]
  (let [;; Generate test data in clinical mass range (2000-20000 Da)
        test-spectra (repeatedly n-spectra
                                 (fn []
                                   (let [start-mz (+ 2000 (rand-int 5000)) ; Random start between 2000-7000
                                         mz-values (range start-mz (+ start-mz n-points) 2) ; 2 Da spacing
                                         ;; Add some realistic peaks with noise
                                         intensities (map-indexed
                                                      (fn [i _]
                                                        (let [base-noise (+ 50 (rand 100))
                                                              ;; Add occasional peaks
                                                              peak (if (and (< 20 i (- n-points 20))
                                                                            (< (rand) 0.05)) ; 5% chance of peak
                                                                     (* 10 (+ 200 (rand 800)))
                                                                     0)]
                                                          (+ base-noise peak)))
                                                      mz-values)]
                                     (create-spectrum
                                      (str "test-" (rand-int 10000))
                                      mz-values
                                      intensities
                                      {:test-data true}))))

        ;; Benchmark clinical preprocessing
        start-time (System/nanoTime)
        clinical-results (mapv clinical-preprocess test-spectra)
        clinical-time (/ (- (System/nanoTime) start-time) 1000000.0)

        ;; Filter out empty spectra (shouldn't happen with realistic data)
        valid-results (filterv #(not (preprocessing/is-empty-spectrum? %)) clinical-results)

        ;; Benchmark peak detection on valid spectra
        start-time2 (System/nanoTime)
        peak-results (mapv #(detect-peaks % :preprocess? false) valid-results)
        peak-time (/ (- (System/nanoTime) start-time2) 1000000.0)

        peak-counts (mapv count peak-results)
        total-peaks (reduce + peak-counts)]

    {:n-spectra n-spectra
     :n-points n-points
     :valid-spectra-after-processing (count valid-results)
     :clinical-preprocessing {:total-time-ms clinical-time
                              :time-per-spectrum-ms (/ clinical-time n-spectra)}
     :peak-detection {:total-time-ms peak-time
                      :time-per-spectrum-ms (if (pos? (count valid-results))
                                              (/ peak-time (count valid-results))
                                              0)}
     :total-peaks-found total-peaks
     :average-peaks-per-spectrum (if (pos? (count valid-results))
                                   (/ total-peaks (count valid-results))
                                   0)
     :peak-counts peak-counts}))

(defn version-info
  "Get MaldiCloj version and system information."
  []
  {:maldicloj-version "0.1.0-SNAPSHOT"
   :clojure-version (clojure-version)
   :java-version (System/getProperty "java.version")
   :available-memory-mb (/ (.maxMemory (Runtime/getRuntime)) 1024 1024)
   :processors (.availableProcessors (Runtime/getRuntime))})

(defn help
  "Display comprehensive help for MaldiCloj core functions."
  []
  (println "
=== MaldiCloj Core API Help ===

🧬 Spectrum Creation (Auto-Optimized):
  (create-spectrum id mz-values intensities metadata)  - Create optimized spectrum
  (spectrum-info spectrum)                             - Get spectrum details + performance info
  (valid-spectrum? spectrum)                           - Validate spectrum

📊 Preprocessing Pipelines:
  (clinical-preprocess spectrum)                       - Weis et al. clinical pipeline
  (basic-preprocess spectrum)                          - General preprocessing
  
🔍 Analysis:
  (detect-peaks spectrum)                              - Peak detection with SNR filtering
  (peak-statistics peaks)                              - Calculate peak statistics
  (compare-spectra spec1 spec2)                        - Compare two spectra
  
⚡ Batch Processing (Auto-Optimized):
  (batch-preprocess spectra)                           - Process multiple spectra efficiently
  (batch-detect-peaks spectra)                         - Peak detection on batches
  
📈 Performance & Benchmarking:
  (profile-performance spectrum)                       - Profile processing bottlenecks
  (show-optimization-impact n-points)                  - Compare optimized vs standard
  (run-comprehensive-benchmark)                        - Full performance test suite
  (estimate-memory-usage n-points)                     - Memory planning
  (suggest-jvm-optimizations)                          - JVM tuning recommendations
  
🗃️ Data Management:
  (spectra->dataset spectra)                           - Convert to dataset format
  (export-results data :csv \"output.csv\")            - Export analysis results
  
ℹ️ Utilities:
  (benchmark-processing n-spectra n-points)           - Quick performance test
  (version-info)                                       - System information
  (help)                                               - This help message

🚀 Performance Features (Automatic):
  • Smart memory layout selection (JVM heap vs native buffers)
  • Automatic data type optimization (float32 for speed, float64 for precision)
  • Intelligent batch sizing based on system resources  
  • Built-in performance monitoring and profiling
  • Zero-configuration optimization for clinical workflows

💡 Tips:
  - create-spectrum() automatically optimizes for your dataset size
  - Use profile-performance() to identify processing bottlenecks
  - Large datasets (>20k points) automatically use native buffers
  - All functions include performance metadata for introspection
  
For detailed documentation, see individual function docstrings.
"))

(defn profile-performance
  "Profile spectrum processing to identify bottlenecks and show optimization impact."
  [spectrum & {:keys [detailed?] :or {detailed? false}}]
  (let [start-time (System/nanoTime)

        ;; Time baseline correction
        baseline-start (System/nanoTime)
        baseline-corrected (basic-preprocess spectrum :baseline-correct? true :smooth? false :normalize? false)
        baseline-time (/ (- (System/nanoTime) baseline-start) 1000000.0)

        ;; Time peak detection
        peak-start (System/nanoTime)
        peaks (detect-peaks baseline-corrected :preprocess? false)
        peak-time (/ (- (System/nanoTime) peak-start) 1000000.0)

        ;; Time full preprocessing
        preproc-start (System/nanoTime)
        preprocessed (clinical-preprocess spectrum)
        preproc-time (/ (- (System/nanoTime) preproc-start) 1000000.0)

        total-time (/ (- (System/nanoTime) start-time) 1000000.0)

        info (spectrum-info spectrum)
        perf-info (:performance info)]

    (println (format "\n🔬 Performance Profile: %s" (:id spectrum)))
    (println "=" (apply str (repeat 50 "=")))
    (println (format "Dataset: %d points (%.1f MB estimated)"
                     (:length info)
                     (get-in info [:memory-usage :estimated-memory-mb] 0.0)))

    (when perf-info
      (println (format "Optimization: %s buffers, %s precision"
                       (if (:native-buffers? perf-info) "Native" "JVM")
                       (name (:mz-dtype perf-info)))))

    (println (format "\n⏱️  Timing Results:"))
    (println (format "  Baseline correction: %.1f ms" baseline-time))
    (println (format "  Peak detection: %.1f ms" peak-time))
    (println (format "  Full preprocessing: %.1f ms" preproc-time))
    (println (format "  Total time: %.1f ms" total-time))

    (println (format "\n📊 Analysis:"))
    (println (format "  Peaks found: %d" (count peaks)))
    (println (format "  Processing rate: %.0f points/ms" (/ (:length info) total-time)))

    (when (and perf-info detailed?)
      (println (format "\n🔧 Performance Details:"))
      (println (format "  Container type: %s" (:container-type perf-info)))
      (println (format "  MZ data type: %s" (:mz-dtype perf-info)))
      (println (format "  Intensity data type: %s" (:intensity-dtype perf-info)))
      (println (format "  Auto-optimized: %s" (:auto-optimized? perf-info))))

    {:total-time-ms total-time
     :baseline-time-ms baseline-time
     :peak-detection-time-ms peak-time
     :preprocessing-time-ms preproc-time
     :peaks-found (count peaks)
     :processing-rate-points-per-ms (/ (:length info) total-time)
     :performance-info perf-info}))

(defn show-optimization-impact
  "Compare performance before and after optimizations."
  [n-points]
  (let [mz-values (range 2000.0 (+ 2000.0 n-points) 1.0)
        intensities (repeatedly n-points #(+ 50 (rand 1000)))

        ;; Create spectrum with old approach (force JVM heap)
        old-spectrum (create-spectrum "unoptimized" mz-values intensities {}
                                      :use-native-buffers? false :mz-dtype :float64)

        ;; Create spectrum with new optimized approach  
        new-spectrum (create-spectrum "optimized" mz-values intensities {})

        ;; Profile both
        old-results (profile-performance old-spectrum)
        new-results (profile-performance new-spectrum)

        speedup (/ (:total-time-ms old-results) (:total-time-ms new-results))
        old-memory (get-in (spectrum-info old-spectrum) [:memory-usage :estimated-memory-mb])
        new-memory (get-in (spectrum-info new-spectrum) [:memory-usage :estimated-memory-mb])
        memory-savings (- old-memory new-memory)]

    (println "\n🚀 Optimization Impact Summary")
    (println "=" (apply str (repeat 40 "=")))
    (println (format "Dataset size: %d points" n-points))
    (println (format "Performance improvement: %.2fx faster" speedup))
    (println (format "Memory savings: %.2f MB (%.1f%% reduction)"
                     memory-savings
                     (* 100 (/ memory-savings old-memory))))
    (println (format "Processing rate: %.0f → %.0f points/ms"
                     (:processing-rate-points-per-ms old-results)
                     (:processing-rate-points-per-ms new-results)))

    {:speedup-factor speedup
     :memory-savings-mb memory-savings
     :memory-reduction-percent (* 100 (/ memory-savings old-memory))
     :old-performance old-results
     :new-performance new-results}))

(defn estimate-memory-usage
  "Estimate memory usage for spectrum processing and analysis."
  [n-points & {:keys [use-native? use-float32?]
               :or {use-native? false use-float32? false}}]
  (let [mz-bytes (* n-points (if use-float32? 4 8))
        intensity-bytes (* n-points 4) ; Always float32 for intensities
        baseline-bytes (* n-points 4) ; Baseline estimation working memory
        peak-overhead (* n-points 0.1) ; Estimated peak storage overhead

        total-bytes (+ mz-bytes intensity-bytes baseline-bytes peak-overhead)
        total-mb (/ total-bytes 1024 1024)]

    {:n-points n-points
     :mz-values-mb (/ mz-bytes 1024 1024)
     :intensities-mb (/ intensity-bytes 1024 1024)
     :processing-overhead-mb (/ (+ baseline-bytes peak-overhead) 1024 1024)
     :total-estimated-mb total-mb
     :use-native-buffers? use-native?
     :precision-mode (if use-float32? "high-speed" "high-precision")}))

(defn suggest-jvm-optimizations
  "Provide JVM optimization recommendations for MaldiCloj."
  []
  (let [runtime (Runtime/getRuntime)
        max-memory-gb (/ (.maxMemory runtime) 1024 1024 1024)
        cores (.availableProcessors runtime)]

    {:system-info {:memory-gb (double max-memory-gb)
                   :cores cores}

     :jvm-flags
     {:heap ["-Xmx8g" "-Xms4g"]
      :gc ["-XX:+UseG1GC" "-XX:MaxGCPauseMillis=200"]
      :performance ["-XX:+UseCompressedOops" "-XX:+OptimizeStringConcat"]
      :native-memory ["-XX:MaxDirectMemorySize=2g"]}

     :batch-recommendations
     (cond
       (< max-memory-gb 8) {:batch-size 25 :parallel-threshold 20000}
       (<= 8 max-memory-gb 16) {:batch-size 50 :parallel-threshold 10000}
       (> max-memory-gb 16) {:batch-size 100 :parallel-threshold 5000})

     :tips ["Use create-spectrum auto-optimization by default"
            "Process spectra in batches for memory efficiency"
            "Profile with profile-performance to find bottlenecks"
            "Use estimate-memory-usage for capacity planning"
            "Monitor JVM memory usage during large batch processing"]}))

(defn run-comprehensive-benchmark
  "Run comprehensive performance benchmarks across multiple dataset sizes."
  [& {:keys [test-sizes include-memory-analysis?]
      :or {test-sizes [1000 5000 10000 25000 50000] include-memory-analysis? false}}]

  (println "\n🔬 MaldiCloj Comprehensive Performance Benchmark")
  (println "=" (apply str (repeat 55 "=")))

  (let [results (atom [])]
    (doseq [size test-sizes]
      (println (format "\nBenchmarking %d points..." size))

      (let [mz-values (range 2000.0 (+ 2000.0 size) 1.0)
            intensities (repeatedly size #(+ 50 (rand 1000)))

            ;; Create optimized spectrum
            start-time (System/nanoTime)
            spectrum (create-spectrum "benchmark" mz-values intensities {})
            creation-time (/ (- (System/nanoTime) start-time) 1000000.0)

            ;; Profile full processing pipeline
            profile-results (profile-performance spectrum)

            ;; Memory analysis
            memory-info (when include-memory-analysis?
                          (estimate-memory-usage size
                                                 :use-native? (get-in profile-results [:performance-info :native-buffers?])
                                                 :use-float32? (= :float32 (get-in profile-results [:performance-info :mz-dtype]))))

            result {:size size
                    :creation-time-ms creation-time
                    :total-processing-time-ms (:total-time-ms profile-results)
                    :processing-rate-points-per-ms (:processing-rate-points-per-ms profile-results)
                    :peaks-found (:peaks-found profile-results)
                    :performance-config (:performance-info profile-results)
                    :memory-analysis memory-info}]

        (swap! results conj result)

        ;; Print results
        (println (format "  Creation: %.1f ms" creation-time))
        (println (format "  Processing: %.1f ms" (:total-processing-time-ms result)))
        (println (format "  Rate: %.0f points/ms" (:processing-rate-points-per-ms result)))
        (println (format "  Config: %s buffers, %s precision"
                         (if (get-in result [:performance-config :native-buffers?]) "Native" "JVM")
                         (name (get-in result [:performance-config :mz-dtype]))))

        (when include-memory-analysis?
          (println (format "  Memory: %.2f MB estimated" (:total-estimated-mb memory-info))))))

    (println "\n✅ Benchmark complete!")

    ;; Return comprehensive results
    {:benchmark-results @results
     :system-info (suggest-jvm-optimizations)
     :summary {:fastest-rate (apply max (map :processing-rate-points-per-ms @results))
               :total-datasets-tested (count test-sizes)
               :auto-optimization-used true}}))
