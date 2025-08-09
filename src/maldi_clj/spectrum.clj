(ns maldi-clj.spectrum
  "High-performance spectrum data structures using dtype-next and Malli validation"
  (:require [malli.core :as m]
            [malli.error :as me]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [tech.v3.dataset :as ds]))

;; ==============================================================================
;; HIGH-PERFORMANCE SPECTRUM SCHEMA - Using dtype-next containers
;; ==============================================================================

;; Malli schema for efficient spectrum with dtype-next containers
(def spectrum-schema
  [:map
   [:id string?]
   [:mz-values [:fn {:error/message "Must be a dtype-next numeric container"}
                #(and (dtype/reader? %)
                      (#{:float32 :float64} (dtype/elemwise-datatype %)))]]
   [:intensities [:fn {:error/message "Must be a dtype-next numeric container"}
                  #(and (dtype/reader? %)
                        (#{:float32 :float64} (dtype/elemwise-datatype %)))]]
   [:metadata map?]])

;; Legacy schema for backward compatibility during migration
(def legacy-spectrum-schema
  [:map
   [:id string?]
   [:mz-values [:vector number?]]
   [:intensities [:vector number?]]
   [:metadata map?]])

;; ==============================================================================
;; CORE SPECTRUM FUNCTIONS - Optimized with dtype-next
;; ==============================================================================

(defn valid-spectrum?
  "Check if spectrum conforms to high-performance schema (with fallback to legacy)"
  [spectrum]
  (or (m/validate spectrum-schema spectrum)
      (m/validate legacy-spectrum-schema spectrum)))

(defn explain-spectrum-errors
  "Get detailed validation errors for spectrum"
  [spectrum]
  (when-let [errors (or (m/explain spectrum-schema spectrum)
                        (m/explain legacy-spectrum-schema spectrum))]
    (me/humanize errors)))

(defn create-spectrum
  "Create an efficient spectrum using dtype-next containers
   
   Args:
     id - String identifier for the spectrum
     mz-values - Collection of m/z values (will be converted to :float64 container)  
     intensities - Collection of intensity values (will be converted to :float32 container)
     metadata - Map of additional spectrum information
   
   Options:
     :mz-dtype - Data type for m/z values (:float32 or :float64, default :float64)
     :intensity-dtype - Data type for intensities (:float32 or :float64, default :float32)
     :container-type - Container type (:jvm-heap or :native-heap, default :jvm-heap)
   
   Returns optimized spectrum with dtype-next containers for high performance"
  [id mz-values intensities metadata & {:keys [mz-dtype intensity-dtype container-type]
                                        :or {mz-dtype :float64
                                             intensity-dtype :float32
                                             container-type :jvm-heap}}]
  (let [;; Create efficient dtype-next containers
        mz-container (dtype/make-container container-type mz-dtype mz-values)
        intensity-container (dtype/make-container container-type intensity-dtype intensities)

        ;; Construct spectrum
        spectrum {:id id
                  :mz-values mz-container
                  :intensities intensity-container
                  :metadata metadata}]

    ;; Validate before returning
    (if (valid-spectrum? spectrum)
      spectrum
      (throw (ex-info "Invalid spectrum data"
                      {:errors (explain-spectrum-errors spectrum)
                       :spectrum spectrum})))))

;; ==============================================================================
;; SPECTRUM UTILITIES - High-performance operations
;; ==============================================================================

(defn spectrum-length
  "Get number of data points in spectrum (zero-copy operation)"
  [spectrum]
  (if (dtype/reader? (:mz-values spectrum))
    (dtype/ecount (:mz-values spectrum))
    (count (:mz-values spectrum))))

(defn get-mz-range
  "Get m/z range as [min, max] using efficient reductions"
  [spectrum]
  (let [mz-data (:mz-values spectrum)]
    (when (pos? (if (dtype/reader? mz-data)
                  (dtype/ecount mz-data)
                  (count mz-data)))
      (if (dtype/reader? mz-data)
        [(dfn/reduce-min mz-data)
         (dfn/reduce-max mz-data)]
        [(apply min mz-data)
         (apply max mz-data)]))))

(defn get-intensity-stats
  "Get intensity statistics using efficient dtype-next operations"
  [spectrum]
  (let [intensities (:intensities spectrum)
        n (if (dtype/reader? intensities)
            (dtype/ecount intensities)
            (count intensities))]
    (when (pos? n)
      (if (dtype/reader? intensities)
        {:count n
         :sum (dfn/sum intensities)
         :mean (dfn/mean intensities)
         :min (dfn/reduce-min intensities)
         :max (dfn/reduce-max intensities)}
        {:count n
         :sum (reduce + intensities)
         :mean (/ (reduce + intensities) n)
         :min (apply min intensities)
         :max (apply max intensities)}))))

(defn spectrum->clojure
  "Convert efficient spectrum back to Clojure collections for compatibility
   Use sparingly - prefer working with dtype-next containers directly"
  [spectrum]
  {:id (:id spectrum)
   :mz-values (if (dtype/reader? (:mz-values spectrum))
                (dtype/->vector (:mz-values spectrum))
                (:mz-values spectrum))
   :intensities (if (dtype/reader? (:intensities spectrum))
                  (dtype/->vector (:intensities spectrum))
                  (:intensities spectrum))
   :metadata (:metadata spectrum)})

(defn clojure->spectrum
  "Convert Clojure collection-based spectrum to efficient dtype-next version"
  [clj-spectrum & options]
  (apply create-spectrum
         (:id clj-spectrum)
         (:mz-values clj-spectrum)
         (:intensities clj-spectrum)
         (:metadata clj-spectrum)
         options))

;; ==============================================================================
;; SPECTRUM COLLECTIONS - Using tech.ml.dataset for multi-spectrum analysis
;; ==============================================================================

(defn spectra->dataset
  "Convert collection of spectra to tech.ml.dataset for efficient analysis
   
   Each spectrum becomes multiple rows (one per m/z, intensity pair) with:
   - :spectrum-id - identifier for grouping  
   - :mz - m/z value
   - :intensity - intensity value
   - Additional columns from metadata"
  [spectra]
  (let [rows (mapcat
              (fn [spectrum]
                (let [{:keys [id mz-values intensities metadata]} spectrum
                      n (if (dtype/reader? mz-values)
                          (dtype/ecount mz-values)
                          (count mz-values))]
                  (map (fn [i]
                         (merge {:spectrum-id id
                                 :mz (if (dtype/reader? mz-values)
                                       (mz-values i)
                                       (nth mz-values i))
                                 :intensity (if (dtype/reader? intensities)
                                              (intensities i)
                                              (nth intensities i))}
                                metadata))
                       (range n))))
              spectra)]
    (ds/->>dataset rows)))

(defn dataset->spectra
  "Convert tech.ml.dataset back to individual spectrum objects"
  [dataset]
  (->> (ds/group-by-column dataset :spectrum-id)
       vals
       (map (fn [group]
              (let [first-row (first (ds/rows group))
                    spectrum-id (:spectrum-id first-row)
                    mz-values (ds/column group :mz)
                    intensities (ds/column group :intensity)
                    ;; Extract metadata (excluding spectrum-specific columns)
                    metadata (dissoc first-row :spectrum-id :mz :intensity)]
                {:id spectrum-id
                 :mz-values mz-values ; Already dtype-next containers
                 :intensities intensities
                 :metadata metadata})))))

;; ==============================================================================
;; COMPATIBILITY LAYER - For existing code migration  
;; ==============================================================================

(defn ^:deprecated create-spectrum-legacy
  "Legacy spectrum creation for backward compatibility
   DEPRECATED: Use create-spectrum for optimal performance"
  [id mz-values intensities metadata]
  {:id id
   :mz-values (vec mz-values) ; Force to vector for legacy compatibility
   :intensities (vec intensities)
   :metadata metadata})

;; Type checking utilities
(defn dtype-spectrum?
  "Check if spectrum uses dtype-next containers"
  [spectrum]
  (and (dtype/reader? (:mz-values spectrum))
       (dtype/reader? (:intensities spectrum))))

(defn legacy-spectrum?
  "Check if spectrum uses legacy Clojure vectors"
  [spectrum]
  (and (vector? (:mz-values spectrum))
       (vector? (:intensities spectrum))))