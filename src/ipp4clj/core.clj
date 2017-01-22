(ns ipp4clj.core)

(use '(incanter stats))
(use 'clojure.core.matrix)
(require '[incanter.core :as ic])
(require '[incanter.charts :as ip])
(require '[clojure.data.avl :as avl])
(require '[clojure.core.matrix.linear :as la])
(use 'clojure.core.matrix.operators)
(require '[clojure.math.combinatorics :as combo])










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
(mean-conditional grid y obs-var 0.4 1.3)
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
