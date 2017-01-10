(ns ipp4clj.core-test
  (:require [clojure.test :refer :all]
            [ipp4clj.core :refer :all]))

(deftest a-test
  (testing "FIXME, I fail."
    (is (= 0 1))))

(deftest single-observation
 "Compare logdensity evaluation with univariate normal")

(defn stationary-covariance
 [gp-variance gp-length-scale]
 (let [delay 1
       X (regression gp-length-scale delay)
       Q (innovation gp-variance gp-length-scale delay)
       P (statcov gp-variance gp-length-scale)
       XPXQ (+ (mmul X P (transpose X)) Q)]
   (mget XPXQ 0 0)))

(defn two-observations
 [times gp-variance gp-length-scale]
 "Covariance between second and first observation should be XQ
  compare this with what is being created from matern-cov-matrix"
 (let [delay (- (second times) (first times))
       X (regression gp-length-scale delay)
       P (statcov gp-variance gp-length-scale)
       K (matern-cov-matrix times gp-variance gp-length-scale)]
   [(mget K 0 1)  (mget (mmul X P) 0 0)]))

(defn matern-cov-matrix
  [times gp-variance gp-length-scale]
  (let [p (ic/length times)
        time-pairs (combo/cartesian-product times times)
        delays (-> (map (partial apply -) time-pairs) abs)
        matern-correlation (fn [d] (let [l (/ (sqrt 5) gp-length-scale)
                                         dl (* d l)]
                                     (* (+ 1 dl (/ (square dl) 3))
                                        (exp (negate dl)))))
        covs (* gp-variance (map matern-correlation delays))]
    (reshape covs [p p])))

(defn trusted-log-likelihood
 "Slow O(n^3) implementation"
 [times observations observation-variance gp-variance gp-length-scale]
 (let [p (ic/length times)
       I (identity-matrix p)
       K (matern-cov-matrix times gp-variance gp-length-scale)
       cov-matrix (+ (* observation-variance I) K)
       exponent (* 0.5 (dot observations (la/solve cov-matrix observations)))
       normalizer (+ (* 0.5 (log (det cov-matrix)))
                     (* p 0.5 (log (* 2 Math/PI))))]
   (negate (+ normalizer exponent))))

(defn trusted-mean-conditional
 [times observations observation-variance gp-variance gp-length-scale]
 (let [p (ic/length times)
       I (identity-matrix p)
       K (matern-cov-matrix times gp-variance gp-length-scale)
       cov-matrix (+ (* observation-variance I) K)
       ones (fill (new-vector p) 1.0)
       prec (dot ones (la/solve cov-matrix ones))
       precxmean (dot ones (la/solve cov-matrix observations))]
   [(/ precxmean prec) (/ 1 prec)]))
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
