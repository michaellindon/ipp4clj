(ns ipp4clj.gaussianprocess
  (:require [incanter.core :as ic]
            [clojure.core.matrix :refer :all]
            [clojure.math.combinatorics :as combo]))

(defn matern-cov-matrix
  [times gp-variance gp-length-scale]
  (let [p (ic/length times)
        time-pairs (combo/cartesian-product times times)
        delays (-> (map (partial apply -) time-pairs) abs)
        matern-correlation (fn [d] (let [l (/ (sqrt 5) gp-length-scale)
                                         dl (* d l)]
                                     (* (+ 1 dl (/ (square dl) 3))
                                        (exp (negate dl)))))
        covs (mul gp-variance (map matern-correlation delays))]
    (reshape covs [p p])))

(defn trusted-log-likelihood
 "Slow O(n^3) implementation"
 [times observations observation-variance gp-variance gp-length-scale]
 (let [p (ic/length times)
       I (identity-matrix p)
       K (matern-cov-matrix times gp-variance gp-length-scale)
       cov-matrix (add (mul observation-variance I) K)
       exponent (* 0.5 (dot observations (mmul ( inverse cov-matrix) observations)))
       normalizer (+ (* 0.5 (log (det cov-matrix)))
                     (* p 0.5 (log (* 2 Math/PI))))]
   (negate (+ normalizer exponent))))

(defn trusted-mean-conditional
 [times observations observation-variance gp-variance gp-length-scale]
 (let [p (ic/length times)
       I (identity-matrix p)
       K (matern-cov-matrix times gp-variance gp-length-scale)
       cov-matrix (add (mul observation-variance I) K)
       ones (fill (new-vector p) 1.0)
       prec (dot ones (mmul (inverse cov-matrix) ones))
       precxmean (dot ones (mmul (inverse cov-matrix) observations))]
   [(/ precxmean prec) (/ 1 prec)]))


;(defn stationary-covariance
; [gp-variance gp-length-scale]
; (let [delay 1
;       X (regression gp-length-scale delay)
;       Q (innovation gp-variance gp-length-scale delay)
;       P (statcov gp-variance gp-length-scale)
;       XPXQ (+ (mmul X P (transpose X)) Q)]
;   (mget XPXQ 0 0)))

;(deftest stationary-covariance
;  (testing "P[0,0] should be gp-variance"
;    (let [gp-var 0.5
;          gp-time-scale 1.5
;          P (statcov gp-var gp-time-scale)
;          P00 (mget P 0 0)
;          difference (- gp-var P00)
;          eps 0.000001]
;      (is (< (abs difference) eps)))))

;(deftest two-observations-covariance
;  (testing "Covariance between second and first observation, XP[0,0], from the
;   state space construction should match K[0,1] from direct Matern construction"
;    (let [times [1.23 1.56]
;          gp-var 0.8
;          gp-time-scale 1.4
;          delay (- (second times) (first times))
;          X (regression gp-time-scale delay)
;          P (statcov gp-var gp-time-scale)
;          K (matern-cov-matrix times gp-var gp-time-scale)
;          K01 (mget K 0 1)
;          XP00 (mget (mmul X P) 0 0)
;          difference (- K01 XP00)
;          eps 0.000001]
;      (is (< (abs difference) eps)))))

;(deftest single-observation-loglikelihood
;  (testing "log-likelihood evaluation for single observation should match
;            direct comparison with univariate normal(0, obs-var+gp-var)"
;    (let [obs [0.3]
;          times [0.6]
;          obs-var 1.5
;          gp-var 0.2
;          gp-time-scale 0.1
;          loglikAR (log-likelihood times obs obs-var gp-var gp-time-scale)
;          loglik (log (pdf-normal (first obs)
;                                  :mean 0
;                                  :sd (sqrt (+ obs-var gp-var))))
;          difference (- loglikAR loglik)
;          eps 0.000001]
;      (is (< (abs difference) eps)))))
;
;(deftest multi-observation-loglikelihood
;  (testing "log-likelihood from state-space against direct construction"
;    (let [times (sort (sample-uniform 10))
;          obs (sample-uniform 10)
;          obs-var (rand)
;          gp-var (rand)
;          gp-t-scale (rand)
;          loglikAR (log-likelihood times obs obs-var gp-var gp-t-scale)
;          loglik (trusted-log-likelihood times obs obs-var gp-var gp-t-scale)
;          difference (- loglikAR loglik)
;          eps 0.000001]
;      (is (< (abs difference) eps)))))
;
; (defn sample-gp
;  "Gaussian Process Function"
;  [gp-variance covariance-kernel]
;  (let [cond-set (atom {})]
;    (fn [t]
;      (let [{times :times Koo :Koo Fo :Fo} @cond-set]
;        (if (empty? times)
;            (let [F-t (sample-normal 1 :sd (sqrt gp-variance))]
;              (do (swap! cond-set assoc :times [t]
;                                        :Koo (matrix [[1]])
;                                        :Fo [F-t])
;                  F-t))
;            3)))))
