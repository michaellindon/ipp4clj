(require '[clojure.test.check :as tc])
(require '[clojure.test.check.generators :as gen])
(require '[clojure.test.check.properties :as prop])


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

(def log-likelihood-test
  (prop/for-all
     [[times obs] (gen/let [n (gen/resize 10 gen/s-pos-int)
                            times (gen/fmap (comp sort distinct) (gen/vector (gen/double* {:infinite? false :NaN? false :min -10.0 :max 10.0}) n))
                            obs (gen/vector (gen/double* {:infinite? false :NaN? false :min -10.0 :max 10.0}) (ic/length times))]
                    [times obs])
      obs-var (gen/double* {:infinite? false :NaN? false :min 0.01 :max 100})
      gp-var (gen/double* {:infinite? false :NaN? false :min 0.01 :max 100})
      gp-time-scale (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})]
   (let [loglik1 (trusted-log-likelihood times obs obs-var gp-var gp-time-scale)
         loglik2 (log-likelihood times obs obs-var gp-var gp-time-scale)
         difference (- loglik1 loglik2)
         eps 1E-10]
     (< (abs difference) eps))))

(def mean-likelihood-test
  (prop/for-all
     [[times obs] (gen/let [n (gen/resize 10 gen/s-pos-int)
                            times (gen/fmap (comp sort distinct) (gen/vector (gen/double* {:infinite? false :NaN? false :min -10.0 :max 10.0}) n))
                            obs (gen/vector (gen/double* {:infinite? false :NaN? false :min -10.0 :max 10.0}) (ic/length times))]
                    [times obs])
      obs-var (gen/double* {:infinite? false :NaN? false :min 0.01 :max 100})
      gp-var (gen/double* {:infinite? false :NaN? false :min 0.01 :max 100})
      gp-time-scale (gen/double* {:infinite? false :NaN? false :min 0.0001 :max 100})]
   (let [[m1 v1] (trusted-mean-conditional times obs obs-var gp-var gp-time-scale)
         [m2 v2] (sample-gp-mean times obs obs-var gp-var gp-time-scale)
         difference1 (- m1 m2)
         difference2 (- v1 v2)
         difference (max difference1 difference2)
         eps 1E-10]
     (< (abs difference) eps))))

(def foo (reductions + 0 (range 0 100000000)))
(def bar (map (partial * 2) foo))

(defn pdf-nst [x mu sigma df] (/ (pdf-t (/ (- x mu) sigma) :df df) sigma))


(gen/sample (gen/set (gen/elements #{1 2 3 4})))
(gen/sample (gen/set (gen/resize 10 gen/pos-int)))

(gen/sample (gen/bind (gen/not-empty (gen/set (gen/resize 10 gen/pos-int)))
                      #(gen/tuple (gen/return %) (gen/set (gen/elements %)))))


(defn logpdf-t [x mu sigma df]
 (let [z (/ (- x mu) sigma)]
  (negate (+ (log sigma) (* (/ (inc df) 2) (log (inc (/ (square z) df))))))))

(defn sample-nst [n mu sigma df] (+ mu (* sigma (sample-t n :df df))))
(defn mean-param [x mu g] (dot (select x g) (select mu g)))


(defn KL [X mu scaler df g]
 (let [full-mean (dot X mu)
       full-sigma (* (square scaler) (inc (dot X X)))
       sel-mean (dot (select X g) (select mu g))
       sel-sigma (* (square scaler) (inc (dot (select X g) (select X g))))
       draws (sample-nst 10000 full-mean (sqrt full-sigma) df)
       full-logpdf (fn [x] (logpdf-t x full-mean (sqrt full-sigma) df))
       sel-logpdf (fn [x] (logpdf-t x sel-mean (sqrt sel-sigma) df))
       difference (fn [x] (- (full-logpdf x) (sel-logpdf x)))]
   (mean (map difference draws))))


(defn KL [X mu scaler g]
 (let [full-mean (dot X mu)
       full-var (* (square scaler) (inc (dot X X)))
       sel-mean (dot (select X g) (select mu g))
       sel-var (* (square scaler) (inc (dot (select X g) (select X g))))
       draws (sample-normal 100000 :mean full-mean :sd (sqrt full-var))
       full-logpdf (fn [x] (logpdf-normal x full-mean full-var))
       sel-logpdf (fn [x] (logpdf-normal x sel-mean sel-var))
       difference (fn [x] (- (full-logpdf x) (sel-logpdf x)))]
   (mean (map difference draws))))
