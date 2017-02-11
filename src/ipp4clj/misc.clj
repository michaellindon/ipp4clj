(ns ipp4clj.misc
  (:require [incanter.stats :refer :all]
            [clojure.core.matrix :refer :all]
            [incanter.core :as ic]
            [incanter.charts :as ip]
            [clojure.core.matrix.linear :as la]
            [clojure.data.avl :as avl]))

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
