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
  "Create a MALDI-TOF spectrum with validation.
   
   Main entry point for spectrum creation with comprehensive validation
   and automatic optimization for high-performance processing."
  [id mz-values intensities & [metadata & options]]
  (let [metadata (or metadata {})
        {:keys [validate? optimize? mz-dtype intensity-dtype]
         :or {validate? true optimize? true
              mz-dtype :float64 intensity-dtype :float32}} options
        spectrum (if optimize?
                   (spectrum/create-spectrum id mz-values intensities metadata
                                             :mz-dtype mz-dtype :intensity-dtype intensity-dtype)
                   (spectrum/create-spectrum-legacy id mz-values intensities metadata))]
    (if validate?
      (if (spectrum/valid-spectrum? spectrum)
        spectrum
        (throw (ex-info "Invalid spectrum created"
                        {:errors (spectrum/explain-spectrum-errors spectrum)
                         :spectrum spectrum})))
      spectrum)))

(defn valid-spectrum?
  "Check if spectrum is valid according to MaldiCloj standards."
  [spectrum]
  (spectrum/valid-spectrum? spectrum))

(defn spectrum-info
  "Get comprehensive information about a spectrum."
  [spectrum]
  {:id (:id spectrum)
   :metadata (:metadata spectrum)
   :length (spectrum/spectrum-length spectrum)
   :mz-range (spectrum/get-mz-range spectrum)
   :intensity-stats (spectrum/get-intensity-stats spectrum)
   :total-ion-current (preprocessing/total-ion-current spectrum)
   :is-empty? (preprocessing/is-empty-spectrum? spectrum)
   :dtype-optimized? (spectrum/dtype-spectrum? spectrum)})

;; =============================================================================
;; High-Level Preprocessing Functions
;; =============================================================================

(defn clinical-preprocess
  "Apply standard clinical preprocessing pipeline (Weis et al. 2020)."
  [spectrum & {:keys [trim-range savgol-window savgol-order snip-iterations tic-normalize?]
               :or {trim-range [2000 20000] savgol-window 21 savgol-order 3
                    snip-iterations 20 tic-normalize? true}}]
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
    (case smooth-method
      :savitzky-golay (smoothing/savitzky-golay 11 3)
      :moving-average (preprocessing/preprocess-spectrum-enhanced
                       :smooth-method :moving-average)
      spectrum)

    baseline-correct?
    (baseline/remove-baseline baseline-method)

    normalize?
    (preprocessing/preprocess-spectrum-enhanced
     :calibrate-method normalize-method)))

;; =============================================================================
;; Peak Detection and Analysis
;; =============================================================================

(defn detect-peaks
  "Detect peaks in a spectrum with sensible defaults."
  [spectrum & {:keys [snr-threshold window-size noise-method preprocess?]
               :or {snr-threshold 3.0 window-size 10 noise-method :mad preprocess? true}}]
  (let [processed-spectrum (if preprocess?
                             (basic-preprocess spectrum)
                             spectrum)]
    (peaks/detect-peaks processed-spectrum
                        :snr-threshold snr-threshold
                        :window-half-size window-size
                        :noise-method noise-method)))

(defn peak-statistics
  "Get comprehensive peak statistics."
  [peaks]
  (peaks/peak-statistics peaks))

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
  "Preprocess multiple spectra with consistent parameters."
  [spectra & {:keys [method] :or {method :basic} :as options}]
  (let [process-fn (case method
                     :clinical #(apply clinical-preprocess % (mapcat identity options))
                     :basic #(apply basic-preprocess % (mapcat identity options)))]
    (process-batch spectra process-fn)))

(defn batch-detect-peaks
  "Detect peaks in multiple spectra."
  [spectra & options]
  (let [detect-fn #(apply detect-peaks % (mapcat identity options))]
    (process-batch spectra detect-fn)))

;; =============================================================================
;; Data Export and Conversion
;; =============================================================================

(defn spectra->dataset
  "Convert spectra to a dataset for analysis."
  [spectra & {:keys [include-peaks? peak-options]
              :or {include-peaks? false peak-options {}}}]
  (let [enhanced-spectra (if include-peaks?
                           (map (fn [spectrum]
                                  (let [peaks (apply detect-peaks spectrum
                                                     (mapcat identity peak-options))]
                                    (assoc spectrum
                                           :peaks peaks
                                           :peak-count (count peaks))))
                                spectra)
                           spectra)]
    (spectrum/spectra->dataset enhanced-spectra)))

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
  (let [test-spectra (repeatedly n-spectra
                                 (fn []
                                   (create-spectrum
                                    (str "test-" (rand-int 10000))
                                    (range 1000 (+ 1000 n-points))
                                    (repeatedly n-points #(rand 1000))
                                    {:test-data true})))

        ;; Benchmark clinical preprocessing
        start-time (System/nanoTime)
        clinical-results (batch-preprocess test-spectra :method :clinical)
        clinical-time (/ (- (System/nanoTime) start-time) 1000000.0)

        ;; Benchmark peak detection
        start-time2 (System/nanoTime)
        peak-results (batch-detect-peaks clinical-results)
        peak-time (/ (- (System/nanoTime) start-time2) 1000000.0)]

    {:n-spectra n-spectra
     :n-points n-points
     :clinical-preprocessing {:total-time-ms clinical-time
                              :time-per-spectrum-ms (/ clinical-time n-spectra)}
     :peak-detection {:total-time-ms peak-time
                      :time-per-spectrum-ms (/ peak-time n-spectra)}
     :peak-counts (map count peak-results)}))

(defn version-info
  "Get MaldiCloj version and system information."
  []
  {:maldicloj-version "0.1.0-SNAPSHOT"
   :clojure-version (clojure-version)
   :java-version (System/getProperty "java.version")
   :available-memory-mb (/ (.maxMemory (Runtime/getRuntime)) 1024 1024)
   :processors (.availableProcessors (Runtime/getRuntime))})

(defn help
  "Display help information for MaldiCloj core functions."
  []
  (println "
=== MaldiCloj Core API Help ===

Main Functions:
  (create-spectrum id mz-values intensities metadata)  - Create a spectrum
  (clinical-preprocess spectrum)                       - Standard clinical preprocessing
  (basic-preprocess spectrum)                          - Basic preprocessing
  (detect-peaks spectrum)                              - Peak detection
  (compare-spectra spec1 spec2)                        - Compare spectra
  
Batch Processing:
  (batch-preprocess spectra)                           - Process multiple spectra
  (batch-detect-peaks spectra)                         - Peak detection on batch
  
Data Management:
  (spectra->dataset spectra)                           - Convert to dataset
  (export-results data :csv \"output.csv\")            - Export results
  
Utilities:
  (spectrum-info spectrum)                             - Get spectrum information
  (benchmark-processing 100 1000)                     - Performance testing
  (version-info)                                       - System information
  (help)                                               - This help message

For detailed documentation, see individual function docstrings.
"))
