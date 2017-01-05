(ns ipp4clj.core)

(use '(incanter stats))
(use 'clojure.core.matrix)
(require '[incanter.core :as ic])
(require '[incanter.charts :as ip])
(require '[clojure.data.avl :as avl])
(require '[clojure.core.matrix.linear :as la])
(use 'clojure.core.matrix.operators)

(defn probit [x] (cdf-normal x))

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




(defn sample-neighbour
    "Sample F[i+1] | F[i] where F[j] is F(t(j)) OR
     Sample F[i-1 | F[i] depending on context]"
    [F regression innovation delay]
    (let [X (regression delay)
          q (sample-safe-mvn (innovation delay))]
     (+ (mmul X F) q)))

(defn sample-sandwich
  "Sample F[i] | F[i-1],F[i+1] where F[j] is F(t(j))"
  [Fv delay-v Fn delay-n regression innovation]
  (let [Xv (regression delay-v)
        Qv (innovation delay-v)
        Xn (regression delay-n)
        delay-nv (+ delay-n delay-v)
        Xnv (regression delay-nv)
        Qnv (innovation delay-nv)
        Qnvi (inverse Qnv)
        XvFv (mmul Xv Fv)
        XnvFv (mmul Xnv Fv)
        QvXn' (mmul Qv (transpose Xn))
        mean-vector (+ XvFv (mmul QvXn' Qnvi (- Fn XnvFv)))
        cov-matrix (- Qv (mmul QvXn' Qnvi (transpose QvXn')))]
    (+ mean-vector (sample-safe-mvn cov-matrix))))

(defn sample-gp
  "Gaussian Process Function"
  [variance length-scale]
  (let [cond-set (atom (avl/sorted-map))
        reg (partial regression length-scale)
        inn (partial innovation variance length-scale)
        rev-reg (partial reverse-regression length-scale)
        rev-inn (partial reverse-innovation variance length-scale)]
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
                       Ft (sample-sandwich Fv delay-v Fn delay-n reg inn)]
                    (do (swap! cond-set assoc t Ft) Ft))
       has-vor (let [[tv Fv] vor
                        delay (- t tv)
                        Ft (sample-neighbour Fv reg inn delay)]
                    (do (swap! cond-set assoc t Ft) Ft))
       has-nach (let [[tn Fn] nach
                        delay (- tn t)
                        Ft (sample-neighbour Fn rev-reg rev-inn delay)]
                    (do (swap! cond-set assoc t Ft) Ft))
       :else (let [Ft (sample-safe-mvn (statcov variance length-scale))]
              (do (swap! cond-set assoc t Ft) Ft))))))))

(defn first-f [f x] (first (f x)))


(defn forward-filter
  "Returns (m+1,M+1) from (m,M)"
  [mM params]
  (let [[m M] mM
        [y X Q] params
        XMXQ (+ Q (mmul X M (transpose X)))
        m-out (+ (mmul X m) (/ (mmul (slice XMXQ 1 0)
                                     (- y (mmul (slice X 0 0) m)))
                               (+ observation-variance (mget XMXQ 0 0))))
        M-out (- XMXQ (/ (outer-product (slice XMXQ 1 0) (slice XMXQ 0 0))
                         (+ observation-variance (mget XMXQ 0 0))))]
      [m-out M-out]))

(defn backward-sample
  "Sample F[i] from F[i+1]"
  [F params]
  (let [[[m M] X Q] params
        XMXQ (+ Q (mmul X M (transpose X)))
        XMXQi (inverse XMXQ)
        MX' (mmul M (transpose X))
        mean-vector (+ m (mmul MX' XMXQi (- F (mmul X m))))
        cov-matrix (- M (mmul MX' XMXQi (transpose MX')))]
      (+ mean-vector (sample-safe-mvn cov-matrix))))

(defn FFBS
  "Forward Filtering Backward Sampling Algorithm"
  [times observations observation-variance variance length-scale]
  (let [delays (map - (rest times) times)
        Qs (map (partial innovation variance length-scale) delays)
        Xs (map (partial regression length-scale) delays)
        minit (let [P (statcov 1 1)] (/ (mmul (slice P 1 0) (first y)) (+ 1 (mget P 0 0))))
        Minit (let [P (statcov 1 1)] (- P
                                        (/ (outer-product (slice P 1 0)
                                                          (slice P 0 0))
                                           (+ 1 (mget P 0 0)))))
        mMs (reductions forward-filter [minit Minit] (map vector (rest y) Xs Qs))]))

(def grid (range 0 100 1))
(def gp (sample-gp 1 1))
(def F (partial first-f gp))
(def y (map (fn [x] (+ (F x)  (sample-normal 1))) grid))
(def delays (map (fn [x] (let [[tn tv] x] (- tn tv))) (map vector (rest grid) grid)))
(def Qs (map (partial innovation 1 1) delays))
(def Xs (map (partial regression 1) delays))
(def minit (let [P (statcov 1 1)] (/ (mmul (slice P 1 0) (first y)) (+ 1 (mget P 0 0)))))
(def Minit (let [P (statcov 1 1)] (- P
                                     (/ (outer-product (slice P 1 0)
                                                       (slice P 0 0))
                                        (+ 1 (mget P 0 0))))))

(def mMs (reductions forward-filter [minit Minit] (map vector (rest y) Xs Qs)))


(def Fn (+ (first (last mMs)) (sample-safe-mvn (second (last mMs)))))
(reverse (reductions backward-sample Fn (reverse (map vector (butlast mMs) Xs Qs))))

(defn countf []
  (let [counter (atom 0)]
    (fn [x] (do (swap! counter inc) @counter))))
