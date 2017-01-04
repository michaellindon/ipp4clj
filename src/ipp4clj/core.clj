(ns ipp4clj.core)

(use '(incanter stats))
(use 'clojure.core.matrix)
(require '[incanter.core :as ic])
(require '[incanter.charts :as ip])
(require '[clojure.data.avl :as avl])
(require '[clojure.core.matrix.linear :as la])

(defn probit [x] (cdFnormal x))

(defn sample-ppp
  "Generates a realization from a 1 dimensional Poisson point process with given
  intensity function and max-intensity over the interval [a,b] by thinning"
  [intensity max-intensity a b]
  (let [n (sample-poisson 1 :lambda (* max-intensity (- b a)))
        times (sample-uniform n :min a :max b)
        acceptance-probability (fn [x] (/ (intensity x) max-intensity))
        acceptance-probabilities (map acceptance-probability times)
        uniforms (sample-uniform n)
        triples (map vector times acceptance-probabilities uniforms)
        reducer (fn [y x]
                 (let [[t p u] x]
                   (if (< u p) (conj y t) y)))]
   (reduce reducer [] triples)))


(defn sample-safe-mvn
  "Safe version of sample-mvn"
  [sigma]
  (let [ncols (column-count sigma)
        z (sample-normal ncols)
        d (la/svd sigma)
        {U :U S :S V* :V*} (la/svd sigma)
        D (diagonal-matrix (map (fn [x] (sqrt (max x 0))) S))]
      (mmul U D z)))

(ic/view (ip/histogram (sample-ppp (comp (partial * 10000) probit sin) 10000 0 100) :nbins 100))
(mmul (:U foo) (diagonal-matrix (:S foo)) (:V* foo))



(defn sample-forward
    "Sample F[i+1] | F[i] where F[j] is F(t(j))"
    [F regression innovation delay]
    (let [X (regression delay)
          q (sample-safe-mvn (innovation delay))]
     (add (mmul X F) q)))

(defn sample-sandwich
  "Sample F[i] | F[i-1],F[i+1] where F[j] is F(t(j))"
  [Fv delay-v Fn delay-n regression innovation]
  (let [Xv (regression delay-v)
        Qv (innovation delay-v)
        Xn (regression delay-n)
        delay-nv (+ delay-n delay-v)
        Xnv (regression delay-nv)
        Qnv (innovation delay-nv)
        XvFv (mmul Xv Fv)
        XnvFv (mmul Xnv Fv)
        QvXn' (mmul Qv (transpose Xn))
        mean-vector (+ XvFv (mmul QvXn' (solve Qnv (- Fn XnvFv))))
        cov-matrix (- Qv (mmul QvXn' (solve Qnv (transpose QvXn'))))]
    (+ mean-vector (sample-safe-mvn cov-matrix))))

(defn sample-gp
  "Gaussian Process Function"
  [variance length-scale]
  (let [cond-set (atom (avl/sorted-map))
        reg (partial regression length-scale)
        inn (partial innovation variance length-scale)]
   (fn [t]
    (if (contains? @cond-set t) (@cond-set t)
     (let [vor (avl/nearest @cond-set < t)
           nach  (avl/nearest @cond-set > t)
           has-vor (not (nil? vor))
           has-nach  (not (nil? nach))
           has-both  (and has-vor has-nach)]
      (cond
       has-both (let [[tv Fv] vor
                       [tn Fn] nach
                       delay-v (- t tv)
                       delay-n (- tn t)
                       Ft (sample-sandwich Fv delay-v Fn delay-n regression innovation)])
       has-vor (let [[tv Fv] vor
                        delay (- t tv)
                        Ft (sample-forward Fv reg inn delay)]
                    (do (swap! cond-set assoc t Ft) Ft))
       has-nach (do 5.0)
       :else (let [Ft (sample-safe-mvn (statcov variance length-scale))]
              (do (swap! cond-set assoc t Ft) Ft))))))))

(defn first-f [f x] (first (f x)))


(defn countf []
  (let [counter (atom 0)]
    (fn [x] (do (swap! counter inc) @counter))))
