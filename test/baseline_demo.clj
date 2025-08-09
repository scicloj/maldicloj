;; Simple demo script
(ns baseline-demo
  (:require [maldi-clj.baseline :as baseline]
            [maldi-clj.baseline-test :as bt]
            [maldi-clj.spectrum :as spectrum]
            [tech.v3.datatype.functional :as dfn]))

(println "=== MALDIcloj Baseline Correction Demo ===")

(let [synthetic (bt/create-synthetic-spectrum-with-baseline 100)
      spectrum (:spectrum synthetic)
      corrected (baseline/remove-baseline spectrum :snip :iterations 25)
      baseline-est (baseline/estimate-baseline spectrum :snip :iterations 25)
      stats (baseline/baseline-correction-stats spectrum corrected baseline-est)]
  (println "Original TIC:" (format "%.2f" (:original-tic stats)))
  (println "Corrected TIC:" (format "%.2f" (:corrected-tic stats)))
  (println "Baseline fraction:" (format "%.3f" (:baseline-fraction stats)))
  (println "Signal retained:" (format "%.3f" (:signal-retained stats)))
  (println "Demo completed successfully!"))
