(ns maldi-clj.smoothing
  "Spectral smoothing functions using Fastmath v3 for high-performance signal processing.
   
   This namespace provides smoothing algorithms essential for MALDI-TOF preprocessing,
   particularly Savitzky-Golay filtering which preserves peak shapes while reducing noise."
  (:require [fastmath.signal :as signal]
            [fastmath.core :as m]
            [fastmath.kernel :as kernel]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [maldi-clj.spectrum :as spectrum]))

(defn savitzky-golay
  "Apply Savitzky-Golay smoothing filter using Fastmath v3.
   
   Savitzky-Golay filtering fits local polynomials to data points, preserving
   peak shapes better than simple moving averages. This is critical for 
   MALDI-TOF spectra where peak preservation is essential.
   
   Parameters:
   - spectrum: A spectrum map with :mz and :intensity keys, or intensity tensor
   - window-size: Size of the smoothing window (must be odd, default: 11)
   - poly-order: Order of polynomial fit (default: 3)
   
   Returns:
   - If spectrum map: spectrum with smoothed intensities
   - If tensor: smoothed intensity tensor
   
   Example:
   (savitzky-golay my-spectrum 21 3)  ; 21-point window, cubic polynomial"
  ([spectrum-or-tensor]
   (savitzky-golay spectrum-or-tensor 11 3))
  ([spectrum-or-tensor window-size]
   (savitzky-golay spectrum-or-tensor window-size 3))
  ([spectrum-or-tensor window-size poly-order]
   (cond
     ;; Handle spectrum map
     (map? spectrum-or-tensor)
     (let [intensities (:intensities spectrum-or-tensor)
           smoothed-intensities (savitzky-golay intensities window-size poly-order)]
       (assoc spectrum-or-tensor :intensities smoothed-intensities))

     ;; Handle intensity tensor directly
     :else
     (let [signal-data (dtype/->reader spectrum-or-tensor)
           savgol-filter (signal/savgol-filter window-size poly-order)
           smoothed (savgol-filter signal-data)]
       (dtype/make-container :jvm-heap :float64 smoothed)))))

(defn savitzky-golay-derivative
  "Calculate derivatives using Savitzky-Golay filter.
   
   Computes numerical derivatives while smoothing, useful for peak detection
   and spectral feature analysis.
   
   Parameters:
   - spectrum-or-tensor: Spectrum map or intensity tensor
   - derivative-order: Order of derivative (1 = first derivative, 2 = second, etc.)
   - window-size: Size of smoothing window (default: 11)
   - poly-order: Order of polynomial fit (default: 3)
   
   Returns:
   - Derivative as tensor or in spectrum map format"
  ([spectrum-or-tensor derivative-order]
   (savitzky-golay-derivative spectrum-or-tensor derivative-order 11 3))
  ([spectrum-or-tensor derivative-order window-size poly-order]
   (cond
     ;; Handle spectrum map
     (map? spectrum-or-tensor)
     (let [intensities (:intensities spectrum-or-tensor)
           derivative-intensities (savitzky-golay-derivative intensities derivative-order window-size poly-order)]
       (assoc spectrum-or-tensor :intensities derivative-intensities))

     ;; Handle intensity tensor directly  
     :else
     (let [signal-data (dtype/->reader spectrum-or-tensor)
           deriv-filter (signal/savgol-filter window-size poly-order derivative-order)
           derivative (deriv-filter signal-data)]
       (dtype/make-container :jvm-heap :float64 derivative)))))

(defn moving-average
  "Apply moving average smoothing using Fastmath.
   
   Simple but fast smoothing method. Less sophisticated than Savitzky-Golay
   but useful for quick noise reduction.
   
   Parameters:
   - spectrum-or-tensor: Spectrum map or intensity tensor
   - window-size: Size of averaging window
   
   Returns:
   - Smoothed spectrum or tensor"
  [spectrum-or-tensor window-size]
  (cond
    ;; Handle spectrum map
    (map? spectrum-or-tensor)
    (let [intensities (:intensities spectrum-or-tensor)
          smoothed-intensities (moving-average intensities window-size)]
      (assoc spectrum-or-tensor :intensities smoothed-intensities))

    ;; Handle intensity tensor directly
    :else
    (let [signal-data (dtype/->reader spectrum-or-tensor)
          ma-filter (signal/moving-average-filter window-size)
          smoothed (ma-filter signal-data)]
      (dtype/make-container :jvm-heap :float64 smoothed))))

(defn gaussian-smooth
  "Apply Gaussian kernel smoothing using Fastmath.
   
   Gaussian smoothing provides a smooth, bell-shaped kernel that gradually
   weights nearby points. Good for noise reduction with minimal peak distortion.
   
   Parameters:
   - spectrum-or-tensor: Spectrum map or intensity tensor
   - bandwidth: Standard deviation of Gaussian kernel (controls smoothing strength)
   
   Returns:
   - Smoothed spectrum or tensor"
  [spectrum-or-tensor bandwidth]
  (cond
    ;; Handle spectrum map
    (map? spectrum-or-tensor)
    (let [intensities (:intensities spectrum-or-tensor)
          smoothed-intensities (gaussian-smooth intensities bandwidth)]
      (assoc spectrum-or-tensor :intensities smoothed-intensities))

    ;; Handle intensity tensor directly
    :else
    (let [signal-data (dtype/->reader spectrum-or-tensor)
          ;; Create Gaussian kernel with 6-sigma window
          window-size (int (* 6 bandwidth))
          gaussian-filter (signal/kernel-smoothing-filter
                           (kernel/kernel :gaussian bandwidth)
                           window-size)
          smoothed (gaussian-filter signal-data)]
      (dtype/make-container :jvm-heap :float64 smoothed))))

(defn variance-stabilization
  "Apply square-root variance stabilization transformation.
   
   Square-root transformation stabilizes variance across intensity ranges,
   required for clinical MALDI-TOF ML pipelines (Weis et al. 2020).
   
   Uses Anscombe transform: sqrt(x + 3/8) for Poisson noise.
   
   Parameters:
   - spectrum-or-tensor: Spectrum map or intensity tensor
   - constant: Anscombe constant (default: 0.375 = 3/8)
   
   Returns:
   - Variance-stabilized spectrum or tensor"
  ([spectrum-or-tensor]
   (variance-stabilization spectrum-or-tensor 0.375))
  ([spectrum-or-tensor constant]
   (cond
     ;; Handle spectrum map
     (map? spectrum-or-tensor)
     (let [intensities (:intensities spectrum-or-tensor)
           stabilized-intensities (variance-stabilization intensities constant)]
       (assoc spectrum-or-tensor :intensities stabilized-intensities))

     ;; Handle intensity tensor directly
     :else
     (dfn/sqrt (dfn/+ spectrum-or-tensor constant)))))

(defn clinical-preprocessing
  "Complete clinical preprocessing pipeline combining variance stabilization and smoothing.
   
   Implements the preprocessing steps from Weis et al. (2020) for AMR prediction:
   1. Square-root variance stabilization  
   2. Savitzky-Golay smoothing (half-window-size 10)
   
   Parameters:
   - spectrum: Spectrum map with :mz and :intensity
   - options: Map with optional keys:
     :variance-constant (default: 0.375)
     :savgol-window (default: 21) 
     :savgol-order (default: 3)
   
   Returns:
   - Preprocessed spectrum map"
  ([spectrum]
   (clinical-preprocessing spectrum {}))
  ([spectrum options]
   (let [{:keys [variance-constant savgol-window savgol-order]
          :or {variance-constant 0.375
               savgol-window 21
               savgol-order 3}} options]
     (-> spectrum
         (variance-stabilization variance-constant)
         (savitzky-golay savgol-window savgol-order)))))
