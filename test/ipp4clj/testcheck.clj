(require '[clojure.test.check :as tc])
(require '[clojure.test.check.generators :as gen])
(require '[clojure.test.check.properties :as prop])


(def sort-idempotent-prop
  (prop/for-all [v (gen/vector gen/int)]
    (= (sort v) (sort (sort v)))))


(tc/quick-check 100 sort-idempotent-prop)

(def stationary-covariance
  (prop/for-all [gp-var (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})
                 gp-time-scale (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})]
    (let [P (statcov gp-var gp-time-scale)
          P00 (mget P 0 0)
          difference (- gp-var P00)
          eps 1E-10]
      (< (abs difference) eps))))

(def two-observations-covariance
  (prop/for-all [rawtimes (gen/vector (gen/double* {:infinite? false :NaN? false :min -10.0 :max 10.0}) 2)
                 gp-var (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})
                 gp-time-scale (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})]
    (let [times (sort rawtimes)
          delay (- (second times) (first times))
          X (regression gp-time-scale delay)
          P (statcov gp-var gp-time-scale)
          K (matern-cov-matrix times gp-var gp-time-scale)
          K01 (mget K 0 1)
          XP00 (mget (mmul X P) 0 0)
          difference (- K01 XP00)
          eps 1E-10]
      (< (abs difference) eps))))

(def single-observation-loglikelihood
  (prop/for-all
    [obs (gen/vector (gen/double* {:infinite? false :NaN? false :min -3.0 :max 3.0}) 1)
     times (gen/vector (gen/double* {:infinite? false :NaN? false :min -10.0 :max 10.0}) 1)
     obs-var (gen/double* {:infinite? false :NaN? false :min 0.01 :max 100})
     gp-var (gen/double* {:infinite? false :NaN? false :min 0.01 :max 100})
     gp-time-scale (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})]
   (let [loglikAR (log-likelihood times obs obs-var gp-var gp-time-scale)
         loglik (log (pdf-normal (first obs)
                                 :mean 0
                                 :sd (sqrt (+ obs-var gp-var))))
         difference (- loglikAR loglik)
         eps 1E-10]
     (< (abs difference) eps))))


(def foo (reductions + 0 (range 0 100000000)))
(def bar (map (partial * 2) foo))
