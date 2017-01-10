(ns ipp4clj.core)

(use '(incanter stats))
(use 'clojure.core.matrix)
(require '[incanter.core :as ic])
(require '[incanter.charts :as ip])
(require '[clojure.data.avl :as avl])
(require '[clojure.core.matrix.linear :as la])
(use 'clojure.core.matrix.operators)
(require '[clojure.math.combinatorics :as combo])
(require '[clojure.tools.trace :as db])

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
 [known-values gp-variance gp-length-scale]
 (let [cond-set (atom (into (avl/sorted-map) known-values))
       reg (partial regression gp-length-scale)
       inn (partial innovation gp-variance gp-length-scale)
       rev-reg (partial reverse-regression gp-length-scale)
       rev-inn (partial reverse-innovation gp-variance gp-length-scale)]
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
          :else (let [Ft (sample-safe-mvn (statcov gp-variance gp-length-scale))]
                  (do (swap! cond-set assoc t Ft) Ft))))))))

(defn first-f [f x] (first (f x)))





(defn FFBS
 "Forward Filtering Backward Sampling Algorithm"
 [times observations observation-variance gp-variance gp-length-scale]
 (let [delays (map - (rest times) times)
       Qs (map (partial innovation gp-variance gp-length-scale) delays)
       Xs (map (partial regression gp-length-scale) delays)
       P (statcov gp-variance gp-length-scale)
       mar-obs-var-1 (+ observation-variance (mget P 0 0))
       minit (/ (* (slice P 1 0) (first observations)) mar-obs-var-1)
       Minit (- P (/ (outer-product (slice P 1 0) (slice P 0 0))
                     mar-obs-var-1))

       forward-filter (fn [mM params]
                        (let [[m M] mM
                              [y X Q] params
                              Xm (mmul X m)
                              XMXQ (+ Q (mmul X M (transpose X)))
                              XMXQH (slice XMXQ 0 0)
                              HXMXQH (mget XMXQ 0 0)
                              mar-obs-var (+ observation-variance HXMXQH)
                              m-out (+ Xm
                                       (/ (mmul XMXQH (- y (mget Xm 0)))
                                          mar-obs-var))

                              M-out (- XMXQ
                                       (/ (outer-product XMXQH XMXQH)
                                          mar-obs-var))]
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


(defn log-likelihood
 "Evaluates the log likelihood"
 [times observations observation-variance gp-variance gp-length-scale]
 (let [delays (map - (rest times) times)
       Qs (map (partial innovation gp-variance gp-length-scale) delays)
       Xs (map (partial regression gp-length-scale) delays)
       P (statcov gp-variance gp-length-scale)
       mar-obs-var-1 (+ observation-variance (mget P 0 0))
       minit (/ (* (slice P 1 0) (first observations)) mar-obs-var-1)
       Minit (- P (/ (outer-product (slice P 1 0) (slice P 0 0))
                     mar-obs-var-1))
       sinit (negate (+ (* 0.5 (log (* Math/PI 2 mar-obs-var-1)))
                        (* (/ 0.5 mar-obs-var-1) (square (first observations)))))
       forward-filter (fn [smM params]
                        (let [[partial-sum m M] smM
                              [y X Q] params
                              Xm (mmul X m)
                              Xm0 (mget Xm 0)
                              XMXQ (+ Q (mmul X M (transpose X)))
                              XMXQH (slice XMXQ 0 0)
                              HXMXQH (mget XMXQ 0 0)
                              mar-obs-var (+ observation-variance HXMXQH)
                              m-out (+ Xm
                                       (/ (mmul XMXQH (- y Xm0)) mar-obs-var))

                              M-out (- XMXQ
                                       (/ (outer-product XMXQH XMXQH) mar-obs-var))
                              residual (- y Xm0)
                              exponent (* (/ 0.5 mar-obs-var) (square residual))
                              logdensity (negate (+ (* 0.5 (log (* mar-obs-var
                                                                   Math/PI 2)))
                                                    exponent))
                              s-out (+ partial-sum logdensity)]
                          [s-out m-out M-out]))]
   (first (reduce forward-filter
                   [sinit minit Minit]
                   (map vector (rest observations) Xs Qs)))))


(defn sample-slice
 "Sample a univariate unnormalized log-density g from position x with length L"
 [g w x]
 (let [y (+ (g x) (negate (sample-exp 1)))
       u (sample-beta 1)
       lower-bound (first (filter (fn [x] (< (g x) y))
                                  (iterate  (fn [x] (- x w))
                                            (- x (* u w)))))
       upper-bound (first (filter (fn [x] (< (g x) y))
                                  (iterate (fn [x] (+ x w))
                                           (+ x (* (- 1 u) w)))))]
   (loop [l lower-bound
          u upper-bound]
     (let [z (first (sample-uniform 1 :min l :max u))]
       (cond
         (> (g z) y) z
         (> z x) (recur l z)
         :else (recur z u))))))

(defn sample-gp-mean
 [times observations observation-variance gp-variance gp-length-scale]
 (let [delays (map - (rest times) times)
       Qs (map (partial innovation gp-variance gp-length-scale) delays)
       Xs (map (partial regression gp-length-scale) delays)
       P (statcov gp-variance gp-length-scale)
       mar-obs-var-1 (+ observation-variance (mget P 0 0))
       minit (/ (* (slice P 1 0) (first observations)) mar-obs-var-1)
       Minit (- P (/ (outer-product (slice P 1 0) (slice P 0 0))
                     mar-obs-var-1))
       Cinit (* (slice P 1 0) (/ (first observations) mar-obs-var-1))
       Dinit (negate (/ (slice P 1 0) mar-obs-var-1))
       forward-filter (fn [acc params]
                        (let [[pm prec C D m M] acc
                              [y X Q] params
                              Xm (mmul X m)
                              Xm0 (mget Xm 0)
                              XMXQ (+ Q (mmul X M (transpose X)))
                              XMXQH (slice XMXQ 0 0)
                              HXMXQH (mget XMXQ 0 0)
                              mar-obs-var (+ observation-variance HXMXQH)
                              prec-out (/ (square (inc (mget (mmul X D) 0))) mar-obs-var)
                              pm-out (/ (* (- y (mget (mmul X C) 0)) (inc (mget (mmul X D) 0))) mar-obs-var)
                              F (/ XMXQH mar-obs-var)
                              C-out (+ (mmul X C) (negate (* F (mget (mmul X C) 0))) (* F y))
                              D-out (- (mmul X D) F (* F (mget (mmul X D) 0)))
                              m-out (+ Xm
                                       (/ (mmul XMXQH (- y Xm0)) mar-obs-var))

                              M-out (- XMXQ
                                       (/ (outer-product XMXQH XMXQH) mar-obs-var))]
                          [pm-out prec-out C-out D-out m-out M-out]))
       ffs (reduce forward-filter
                             [(/ (first observations) mar-obs-var-1) (/ 1 mar-obs-var-1) Cinit Dinit minit Minit]
                             (map vector (rest observations) Xs Qs))
       precision (second ffs)
       precxmean (first ffs)]
   [(/ precxmean precision) (/ 1 precision)]))






(def grid (range 1 10 1))
(def gp (sample-gp {} 1 1))
(def F (partial first-f gp))
(def observation-variance 0.3)
(def y (map (fn [x] (+ (F x)  (sample-normal 1 :mean 0 :sd (sqrt observation-variance)))) grid))
(/ (log-likelihood grid y observation-variance 1.4 4.4)
   (log-likelihood grid y observation-variance 10 10))
(/ (trusted-log-likelihood grid y observation-variance 1.4 4.4)
   (trusted-log-likelihood grid y observation-variance 10 10))
(log (pdf-normal (first y) :mean 0 :sd (sqrt (+ observation-variance 10))))
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






(variance (take 10000 (iterate (partial sample-slice (comp log pdf-exp) 1) 0)))
(ic/view (ip/histogram (take 10000 (iterate (partial sample-slice (comp log (fn [x] (pdf-normal x :mean 10 :sd 3))) 0.2) 0.5)) :nbins 100))
(ic/view (ip/histogram (sample-normal 10000 :mean 10 :sd 3) :nbins 100))


(def grid (range 1 10 7))
(def gp (sample-gp {} 1 1))
(def F (partial first-f gp))
(def observation-variance 1.0)
(def y (map (fn [x] (+ (F x) 10  (sample-normal 1 :mean 0 :sd (sqrt observation-variance)))) grid))
(sample-gp-mean grid y observation-variance 1 1)
(trusted-mean-conditional grid y observation-variance 1 1)
