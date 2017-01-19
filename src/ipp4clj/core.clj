(ns ipp4clj.core)

(use '(incanter stats))
(use 'clojure.core.matrix)
(require '[incanter.core :as ic])
(require '[incanter.charts :as ip])
(require '[clojure.data.avl :as avl])
(require '[clojure.core.matrix.linear :as la])
(use 'clojure.core.matrix.operators)
(require '[clojure.math.combinatorics :as combo])


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
 [known-values gp-var gp-time-scale]
 (let [cond-set (atom (into (avl/sorted-map) known-values))
       reg (partial regression gp-time-scale)
       inn (partial innovation gp-var gp-time-scale)
       rev-reg (partial reverse-regression gp-time-scale)
       rev-inn (partial reverse-innovation gp-var gp-time-scale)]
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
          :else (let [Ft (sample-safe-mvn (statcov gp-var gp-time-scale))]
                  (do (swap! cond-set assoc t Ft) Ft))))))))




(defn AR-params
  [times obs obs-var gp-var gp-time-scale]
  (let [delays (map - (rest times) times)
        Qs (map (partial innovation gp-var gp-time-scale) delays)
        Xs (map (partial regression gp-time-scale) delays)
        P (statcov gp-var gp-time-scale)
        mar-obs-var-1 (+ obs-var (mget P 0 0))
        minit (/ (* (slice P 1 0) (first obs)) mar-obs-var-1)
        HP (slice P 0 0)
        Minit (- P (/ (outer-product HP HP) mar-obs-var-1))]
    [delays Qs Xs P HP mar-obs-var-1 minit Minit]))


(defn forward-filter
 [inputs params]
 (let [[m M _ _ _] inputs
       [y obs-var X Q] params
       Xm (mmul X m)
       mar-obs-mean (mget Xm 0)
       XMXQ (+ Q (mmul X M (transpose X)))
       XMXQH (slice XMXQ 0 0)
       HXMXQH (mget XMXQ 0 0)
       mar-obs-var (+ obs-var HXMXQH)
       residual (- y mar-obs-mean)
       m-out (+ Xm (/ (mmul XMXQH residual) mar-obs-var))
       M-out (- XMXQ (/ (outer-product XMXQH XMXQH) mar-obs-var))]
   [m-out M-out mar-obs-mean mar-obs-var XMXQH]))

(defn backward-compose
 [m M X Q]
 (let [XMXQ (+ Q (mmul X M (transpose X)))
       XMXQi (inverse XMXQ)
       MX' (mmul M (transpose X))
       XM (mmul X M)
       Xm (mmul X m)
       mean-function (fn [F] (+ m (mmul MX' XMXQi (- F Xm))))
       cov-matrix (- M (mmul MX' XMXQi XM))]
   [mean-function cov-matrix]))

(defn backward-sample
 [F distribution]
 (let [[mean-function cov-matrix] distribution]
   (+ (mean-function F) (sample-safe-mvn cov-matrix))))

(defn FFBC
 "Forward Filtering Backward Compose Algorithm"
 [times obs obs-var gp-var gp-time-scale]
 (let [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
       (AR-params times obs obs-var gp-var gp-time-scale)
       FF (reductions forward-filter
                      [minit Minit 0 mar-obs-var-1 HP]
                      (map vector (rest obs) (repeat obs-var) Xs Qs))
       BCn [(first (last FF)) (second (last FF))]
       BC (map backward-compose (map first FF) (map second FF) Xs Qs)]
   [BC BCn]))

(defn FFBS
 "Forward Filtering Backward Compose Algorithm"
 [times obs obs-var gp-var gp-time-scale]
 (let [[BC BCn] (FFBC times obs obs-var gp-var gp-time-scale)
       Fn (+ (first BCn) (sample-safe-mvn (second BCn)))
       Fs (reverse (reductions backward-sample Fn (reverse BC)))]
   (zipmap times Fs)))

(defn logpdf-gen-F
 [F times obs obs-var gp-var gp-time-scale]
 (let [[BC BCn] (FFBC times obs obs-var gp-var gp-time-scale)
       Fs (map F times)
       Fn (last Fs)
       mean-functions (map first BC)
       means (map (fn [f x] (f x)) mean-functions (rest Fs))
       covars (map second BC)]
   (+ (reduce + (map logpdf-mvnormal (drop-last Fs) means covars))
      (logpdf-mvnormal Fn (first BCn) (second BCn)))))

(defn logpdf-fc-F
 [F times obs obs-var gp-var gp-time-scale]
 (let [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
       (AR-params times obs obs-var gp-var gp-time-scale)
       Fs (map F times)
       means (cons [0 0 0] (map (fn [X F] (mmul X F)) Xs (drop-last Fs)))
       covars (cons P Qs)
       logdensity-f (reduce + (map logpdf-mvnormal Fs means covars))
       logdensity-y (reduce + (map logpdf-normal obs (map first Fs) (repeat obs-var)))]
   (+ logdensity-y logdensity-f)))


(defn log-likelihood
 "Evaluates the log likelihood"
 [times obs obs-var gp-var gp-time-scale]
 (if (or (neg? gp-var) (neg? gp-time-scale) (neg? obs-var))
     (Double/NEGATIVE_INFINITY)
     (let [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
           (AR-params times obs obs-var gp-var gp-time-scale)
           FF (reductions forward-filter
                          [minit Minit 0 mar-obs-var-1 HP]
                          (map vector (rest obs) (repeat obs-var) Xs Qs))
           mar-means (map (fn [x] (nth x 2)) FF)
           mar-vars (map (fn [x] (nth x 3)) FF)]
       (reduce + (map (fn [x] (let [[y mu s2] x] (logpdf-normal y mu s2)))
                      (map vector obs mar-means mar-vars))))))



(defn sample-gp-mean
 [times obs obs-var gp-var gp-time-scale]
 (let
  [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
   (AR-params times obs obs-var gp-var gp-time-scale)
   Cinit (* (slice P 1 0) (/ (first obs) mar-obs-var-1))
   Dinit (negate (/ (slice P 1 0) mar-obs-var-1))
   forward-filter
    (fn [acc params]
     (let
      [[pm prec C D m M] acc
       [y X Q] params
       [m+1 M+1 _ mar-obs-var XMXQH] (forward-filter [m M 0 0 0] [y obs-var X Q])
       XD (mmul X D)
       HXD (mget XD 0)
       HXD+1 (inc HXD)
       XC (mmul X C)
       HXC (mget XC 0)
       prec+1 (+ prec (/ (square HXD+1) mar-obs-var))
       pm+1 (+ pm (/ (* (- y HXC) HXD+1) mar-obs-var))
       F (/ XMXQH mar-obs-var)
       C+1 (+ XC (negate (* F HXC)) (* F y))
       D+1 (- XD F (* F HXD))]
      [pm+1 prec+1 C+1 D+1 m+1 M+1]))

   ffs (reduce forward-filter
               [(/ (first obs) mar-obs-var-1) (/ 1 mar-obs-var-1) Cinit Dinit minit Minit]
               (map vector (rest obs) Xs Qs))
   precision (second ffs)
   precxmean (first ffs)]
  [(/ precxmean precision) (/ 1.0 precision)]))







(def times (range 1 10 3))
(def gp-var 1)
(def gp-time-scale 1)
(def gp (sample-gp {} gp-var gp-time-scale))
(def F1 (sample-gp {} gp-var gp-time-scale))
(def F2 (sample-gp {} gp-var gp-time-scale))
(def F (comp first gp))
(def obs-var 0.1)
(def obs (map (fn [x] (+ (F x)  (sample-normal 1 :mean 0 :sd (sqrt obs-var)))) times))
(def ARparams (AR-params times obs obs-var gp-var gp-time-scale))
(def delays (first ARparams))
(def Qs (second ARparams))
(def Xs (nth ARparams 2))
(def P (nth ARparams 3))
(def HP (nth ARparams 4))
(def mar-obs-var-1 (nth ARparams 5))
(def minit (nth ARparams 6))
(def Minit (nth ARparams 7))
ARparams [[delays Qs Xs P HP mar-obs-var-1 minit Minit]]

(- (logpdf-fc-F F1 times obs obs-var gp-var gp-time-scale)
   (logpdf-fc-F F2 times obs obs-var gp-var gp-time-scale))
(- (logpdf-gen-F F1 times obs obs-var gp-var gp-time-scale)
   (logpdf-gen-F F2 times obs obs-var gp-var gp-time-scale))

(log-likelihood times obs obs-var gp-var gp-time-scale)
(trusted-log-likelihood times obs obs-var gp-var gp-time-scale)
(log (pdf-normal (first y) :mean 0 :sd (sqrt (+ obs-var 10))))
(def pivots (FFBS times obs obs-var gp-var gp-time-scale))
(def F1 (comp first (sample-gp pivots gp-var gp-time-scale)))
(def F2 (comp first (sample-gp pivots gp-var gp-time-scale)))
(def F3 (comp first (sample-gp pivots gp-var gp-time-scale)))
(def F4 (comp first (sample-gp pivots gp-var gp-time-scale)))
(def F5 (comp first (sample-gp pivots gp-var gp-time-scale)))
(doto (ip/scatter-plot times obs)
  (ip/add-function F1 0 10)
  (ip/add-function F2 0 10)
  (ip/add-function F3 0 10)
  (ip/add-function F4 0 10)
  (ip/add-function F5 0 10)
  ic/view)


(def intensity (comp (partial * max-intensity) probit F))

(defn integrate [f a b]
  (let [points (range a b (/ 10000))]
    (/ (reduce + (map f points)) 10000)))
(riemann-lower-sum F 0 1 10)
(riemann-lower-sum F 0 1 100)
(riemann-lower-sum F 0 1 1000)
(riemann-lower-sum F 0 1 10000)


(ic/view (ip/function-plot intensity 3 4))
(ic/view (ip/function-plot F 3 4))
(def max-intensity 100)
(def intensity (comp (partial * max-intensity) probit F))
(var (take 10000 (iterate (partial sample-slice (comp log pdf-exp) 1) 0)))
(ic/view (ip/histogram (take 10000 (iterate (partial sample-slice (comp log (fn [x] (pdf-normal x :mean 10 :sd 3))) 0.2) 0.5)) :nbins 100))
(ic/view (ip/histogram (sample-normal 10000 :mean 10 :sd 3) :nbins 100))


(def grid (range 1 30 1))
(def gp (sample-gp {} 1 1))
(def F (comp first gp))
(def obs-var 0.5)
(def y (map (fn [x] (+ (F x) -100  (sample-normal 1 :mean 0 :sd (sqrt obs-var)))) grid))
(sample-gp-mean grid y obs-var 0.4 1.3)
(trusted-mean-conditional grid y obs-var 0.4 1.3)


(- (logpdf-fc-F F1 [0.00048828125 0.0009765625 1] [1 1 1] 1 1 1)
   (logpdf-fc-F F2 [0.00048828125 0.0009765625 1] [1 1 1] 1 1 1))
(- (logpdf-gen-F F1 [0.00048828125 0.0009765625 1] [1 1 1] 1 1 1)
   (logpdf-gen-F F2 [0.00048828125 0.0009765625 1] [1 1 1] 1 1 1))


(- (logpdf-fc-F F1 [0.921875 1] [1 1] 1 1 1.5)
   (logpdf-fc-F F2 [0.921875 1] [1 1] 1 1 1.5))
(- (logpdf-fc-F F1 [0.921875 1] [1 1] 1 1 1.5)
   (logpdf-fc-F F2 [0.921875 1] [1 1] 1 1 1.5))

[-0.14739990234375 -0.125 1] [1 1 1] 1 1 8
