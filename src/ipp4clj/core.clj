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
 [known-values variance length-scale]
 (let [cond-set (atom (into (avl/sorted-map) known-values))
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





(defn FFBS
 "Forward Filtering Backward Sampling Algorithm"
 [times observations observation-variance variance length-scale]
 (let [delays (map - (rest times) times)
       Qs (map (partial innovation variance length-scale) delays)
       Xs (map (partial regression length-scale) delays)
       P (statcov 1 1)
       minit (/ (* (slice P 1 0) (first observations)) (+ observation-variance
                                                          (mget P 0 0)))
       Minit (- P (/ (outer-product (slice P 1 0) (slice P 0 0))
                     (+ observation-variance (mget P 0 0))))

       forward-filter (fn [mM params]
                        (let [[m M] mM
                              [y X Q] params
                              Xm (mmul X m)
                              XMXQ (+ Q (mmul X M (transpose X)))
                              XMXQH (slice XMXQ 0 0)
                              HXMXQH (mget XMXQ 0 0)
                              m-out (+ Xm
                                       (/ (mmul XMXQH (- y (mget Xm 0)))
                                          (+ observation-variance HXMXQH)))

                              M-out (- XMXQ
                                       (/ (outer-product XMXQH XMXQH)
                                          (+ observation-variance HXMXQH)))]
                          [m-out M-out]))

       backward-sample (fn [F params]
                           (let [[[m M] X Q] params
                                 XMXQ (+ Q (mmul X M (transpose X)))
                                 XMXQi (inverse XMXQ)
                                 MX' (mmul M (transpose X))
                                 XM (mmul X M)
                                 Xm (mmul X m)
                                 mean-vector (+ m (mmul MX' XMXQi (- F Xm)))
                                 cov-matrix (- M (mmul MX' XMXQi XM))]
                             (+ mean-vector (sample-safe-mvn cov-matrix))))

       mMs (reductions forward-filter
                       [minit Minit]
                       (map vector (rest observations) Xs Qs))
       Fn (+ (first (last mMs)) (sample-safe-mvn (second (last mMs))))
       Fs (reverse (reductions backward-sample
                               Fn (reverse (map vector (butlast mMs) Xs Qs))))]
   (zipmap times Fs)))

(def grid (range 0 10 3))
(def gp (sample-gp {} 1 1))
(def F (partial first-f gp))
(def observation-variance 0.000001)
(def y (map (fn [x] (+ (F x)  (sample-normal 1 :mean 0 :sd (sqrt observation-variance)))) grid))
(def pivots (FFBS grid y observation-variance 1 1))
(def F1 (partial first-f (sample-gp pivots 1 1)))
(def F2 (partial first-f (sample-gp pivots 1 1)))
(def F3 (partial first-f (sample-gp pivots 1 1)))
(def F4 (partial first-f (sample-gp pivots 1 1)))
(def F5 (partial first-f (sample-gp pivots 1 1)))
(doto (ip/scatter-plot grid y)
  (ip/add-function F1 0 10)
  (ip/add-function F2 0 10)
  (ip/add-function F3 0 10)
  (ip/add-function F4 0 10)
  (ip/add-function F5 0 10)
  ic/view)

(def gitter (range 0 10 0.01))
(proto-repl-charts.charts/custom-chart
  "Custom"
  {:data {:columns
          [(cons "x" gitter)
           (cons "F1" (map F1 gitter))
           (cons "F2" (map F2 gitter))
           (cons "F3" (map F3 gitter))
           (cons "F4" (map F4 gitter))
           (cons "F5" (map F5 gitter))]
          :point {:show false}
          :x "x"}})



(def delays (map (fn [x] (let [[tn tv] x] (- tn tv))) (map vector (rest grid) grid)))
(def Qs (map (partial innovation 1 1) delays))
(def Xs (map (partial regression 1) delays))
(def minit (let [P (statcov 1 1)] (/ (mmul (slice P 1 0) (first y)) (+ 1 (mget P 0 0))
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

(ic/view (ip/histogram (sample-ppp (comp (partial * 10000) probit sin) 10000 0 100) :nbins 100))


(proto-repl-charts.charts/custom-chart
  "Custom"
  {:data { :type "line" :columns [["foo" [1 2 3]]]}})

(proto-repl-charts.charts/custom-chart
  "Custom"
  {:data {:columns
          [(cons "x" grid)
           (cons "gp" (map (comp probit F) grid))
           (cons "gp+" (map (fn [x] (probit (+ (F x) 0.5))) grid))
           (cons "gp-" (map (fn [x] (probit (- (F x) 0.5))) grid))]
          :types {:gp "spline" :gp+ "area" :gp- "area"}
          :colors {:gp- "grey" :gp+ "grey"}
          :x "x"}})

(proto-repl-charts.charts/custom-chart
  "Custom"
  {:data {:columns
          [(cons "x" grid)
           (cons "gp" (map (comp probit F) grid))
           (cons "gp+" (map (fn [x] (probit (+ (F x) 0.5))) grid))
           (cons "gp-" (map (fn [x] (probit (- (F x) 0.5))) grid))]
          :types {:gp "spline" :gp+ "area" :gp- "area"}
          :x "x"
          :groups [["gp" "gp+" "gp-"]]}})
