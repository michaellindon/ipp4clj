(defn Q
 "Returns unscaled Covariance matrix for length scale ls and delay d"
  [ls d]
  (let [l (/ (sqrt 5) ls)
        q (/ (* 16 (pow l 5)) 3)
        em2dl (exp (* l d -2))]
    (matrix [
             [(/ (* q (+ 3 (* em2dl (+ -3 (* -2 d l (+ 3 (* d l (+ 3 (* d l (+ 2 (* d l))))))))))) (* 16 (pow l 5)))
              (* 0.125 em2dl q (pow d 4))
              (* -1 (/ (* q (+ (* -1 em2dl) 1 (* 2 d l em2dl (+ -1 (* d l (+ -1 (* d l (+ -2 (* d l))))))))) (* 16 (pow l 3))))]
             [(* 0.125 em2dl q (pow d 4))
              (/ (* q (+ 1 (* em2dl (+ -1 (* -2 d l (+ 1 (* d l (pow (+ -1 (* d l)) 2)))))))) (* 16 (pow l 3)))
              (* 0.125 em2dl q (pow d 2) (pow (+ -2 (* d l)) 2))]
             [(* -1 (/ (* q (+ (* -1 em2dl) 1 (* 2 d l em2dl (+ -1 (* d l (+ -1 (* d l (+ -2 (* d l))))))))) (* 16 (pow l 3))))
              (* 0.125 em2dl q (pow d 2) (pow (+ -2 (* d l)) 2))
              (/ (* q (+ 3 (* em2dl (+ -3 (* -2 d l (+ -5 (* d l (+ 11 (* d l (+ -6 (* d l))))))))))) (* 16 l))]])))

(defn innovation [ρ² l Δ] (mmul ρ² (Q l Δ)))

(defn regression
  "Returns a regression matrix for length scale ls and delay d"
  [ls d]
  (let [l (/ (sqrt 5) ls)
        emdl (exp (* l d -1))]
    (matrix [
             [(* 0.5 emdl (+ 2 (* 2 d l) (* d d l l)))
              (* emdl d (+ 1 (* d l)))
              (* 0.5 emdl d d)]
             [(* -1 0.5 emdl d d l l l)
              (* -1 emdl (+ -1 (* -1 d l) (* d d l l)))
              (* -1 0.5 emdl d (+ -2 (* d l)))]
             [(* 0.5 emdl d l l l (+ -2 (* d l)))
              (* emdl d l l (+ -3 (* d l)))
              (* 0.5 emdl (+ 2 (* -4 d l) (* d d l l)))]])))


(defn statcorr
  "Returns the stationary correlation matrix for length scale ls"
  [ls]
  (let [l (/ (sqrt 5) ls)
        q (/ (* 16 (pow l 5)) 3)]
      (matrix [
               [(/ (* 3 q) (* 16 (pow l 5)))
                0
                (* -1 (/ q (* 16 l l l)))]
               [0
                (/ q (* 16 l l l))
                0]
               [(* -1 (/ q (* 16 l l l)))
                0
                (/ (* 3 q) (* 16 l))]])))

(defn statcov [ρ² l] (mmul ρ² (statcorr l)))
