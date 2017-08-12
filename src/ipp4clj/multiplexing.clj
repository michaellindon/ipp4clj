(ns ipp4clj.multiplexing
  (:require [clojure.core.matrix :refer :all]
            [incanter.stats :refer :all]
            [incanter.core :as ic]
            [clojure.data.avl :as avl]
            [ssm4clj.core :refer :all]
            [ssm4clj.misc :refer :all]
            [ipp4clj.misc :refer :all]
            [clojure.walk :refer :all]
            [distributions.core :as d]
            [clojure.data.json :as json]))

;Change this to Logistic for Polya-Gamma
(def link probit)
(def base-measure (d/normal 0 1))
(def mean-prior-var-prior (d/inverse-gamma 2 0.1))
(defn sample-bool [p] (< (rand) p))
(defn sample-bernoulli [p] (if (sample-bool p) 1 0))
(defn zero-fn [x] 0.0)
(defn zero-3fn [x] (matrix [0.0 0.0 0.0]))
(defn group-bool [pred coll]
  (let [bool-map (group-by pred coll)
        trues (get bool-map true)
        falses (get bool-map false)]
    [trues (if (nil? falses) [] falses)]))

(defn randomize-single-trials [C]
  (let [{max-intensity :max-intensity
         trials :trials} C
        gp-mean (sample-normal 1)
        gp-var (rand)
        gp-time-scale (rand)
        F (sample-gp {} gp-var gp-time-scale)
        f (comp first F)
        intensity (comp (partial * max-intensity) link (partial + gp-mean ) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))
        t-mod (fn [t]
                (let [{obs-times :obs-times} t
                      aug-times (sample-ppp thin-intensity max-intensity 0 1)
                      obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
                      aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
                      obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
                      aug-map (into (avl/sorted-map) (map vector aug-times aug-y))
                      ]
                  {:obs-times obs-times
                   :aug-times aug-times
                   :y (merge obs-map aug-map)
                   :start-t 0
                   :end-t 1} ))
        t-new (map t-mod trials)]
    (assoc C :trials t-new :gp-mean gp-mean :gp-var gp-var :gp-time-scale gp-time-scale :F F)))


(defn generate-single-trials [num-trials gp-mean gp-var gp-time-scale max-intensity]
  (let [F (sample-gp {} gp-var gp-time-scale)
        ;F (fn [x] [-2.5 -2.5 -2.5])
        f (comp first F)
        intensity (comp (partial * max-intensity) link (partial + gp-mean ) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))]
    {:F F :gp-mean gp-mean :gp-var gp-var :gp-time-scale gp-time-scale :max-intensity max-intensity
     :trials
     (take
      num-trials
      (repeatedly (fn []
                    (let [obs-times (sample-ppp intensity max-intensity 0 1)
                          aug-times (sample-ppp thin-intensity max-intensity 0 1)
                          ;obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
                          ;aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
                          obs-y (map (fn [t] (let [w (d/sample (d/polya-gamma (+ gp-mean (f t))))] [(/ 0.5 w) (/ 1 w)])) obs-times)
                          aug-y (map (fn [t] (let [w (d/sample (d/polya-gamma (+ gp-mean (f t))))] [(/ -0.5 w) (/ 1 w)])) aug-times)
                          obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
                          aug-map (into (avl/sorted-map) (map vector aug-times aug-y))]
                        {:obs-times (sort obs-times)
                         :aug-times aug-times
                         :y (merge obs-map aug-map)
                         :start-t 0
                         :end-t 1}))))}))

(defn intensity-functions
  "Provides all relevant intensity functions for AB trials.
  A modular function useful for unit testing"
  [alpha intensity-A intensity-B max-intensity-A max-intensity-B]
  (let [i-Aar (fn [x] (* (- 1 (alpha x)) (intensity-A x)))
        i-Arr (fn [x] (- max-intensity-A (intensity-A x)))
        i-Bar (fn [x] (* (alpha x) (intensity-B x)))
        i-Brr (fn [x] (- max-intensity-B (intensity-B x )))
        i-Aaa (fn [x] (* (alpha x) (intensity-A x)))
        i-Baa (fn [x] (* (- 1 (alpha x)) (intensity-B x)))]
    [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr]))

(defn randomize-dual-trials [A B AB]
  (let [
        {F-A :F gp-mean-A :gp-mean max-intensity-A :max-intensity} A
        {F-B :F gp-mean-B :gp-mean max-intensity-B :max-intensity} B
        f-A (comp first F-A)
        m-A (comp (partial + gp-mean-A) f-A)
        intensity-A (comp (partial * max-intensity-A) link m-A)
        f-B (comp first F-B)
        m-B (comp (partial + gp-mean-B ) f-B)
        intensity-B (comp (partial * max-intensity-B) link m-B)
        gp-time-scale (rand)
        gp-var (rand)
        {trials :trials} AB
        t-mod (fn [t]
                (let [
                      switching (sample-bernoulli 0.5)
                      G (if (= switching 1) (sample-gp {} gp-var gp-time-scale) zero-3fn)
                      g (comp first G)
                      gp-mean-AB (sample-normal 1)
                      alpha (comp link (partial + gp-mean-AB) g)
                      [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr] (intensity-functions
                                                             alpha
                                                             intensity-A
                                                             intensity-B
                                                             max-intensity-A
                                                             max-intensity-B)
                      start-t 0
                      end-t 1
                      Aaa (sample-ppp i-Aaa max-intensity-A start-t end-t)
                      Aar (sample-ppp i-Aar max-intensity-A start-t end-t)
                      Arr (sample-ppp i-Arr max-intensity-A start-t end-t)
                      Baa (sample-ppp i-Baa max-intensity-B start-t end-t)
                      Bar (sample-ppp i-Bar max-intensity-B start-t end-t)
                      Brr (sample-ppp i-Brr max-intensity-B start-t end-t)
                      obs-times (sort (concat Aaa Baa))
                      m-AB (comp (partial + gp-mean-AB) g)
                      Aar-y-A (map (fn [t] (sample-right-normal (m-A t) 1 )) Aar)
                      Aar-y-AB (map (fn [t] (sample-left-normal (m-AB t) 1 )) Aar)
                      Arr-y-A (map (fn [t] (sample-left-normal (m-A t) 1 )) Arr)
                      Bar-y-B (map (fn [t] (sample-right-normal (m-B t) 1 )) Bar)
                      Bar-y-AB (map (fn [t] (sample-right-normal (m-AB t) 1 )) Bar)
                      Brr-y-B (map (fn [t] (sample-left-normal (m-B t) 1 )) Brr)
                      Aaa-y-A (map (fn [t] (sample-right-normal (m-A t) 1 )) Aaa)
                      Aaa-y-AB (map (fn [t] (sample-right-normal (m-AB t) 1 )) Aaa)
                      Baa-y-B (map (fn [t] (sample-right-normal (m-B t) 1 )) Baa)
                      Baa-y-AB (map (fn [t] (sample-left-normal (m-AB t) 1 )) Baa)
                      y-A (into (avl/sorted-map) (concat (map vector Aar Aar-y-A)
                                                         (map vector Arr Arr-y-A)
                                                         (map vector Aaa Aaa-y-A)))
                      y-B (into (avl/sorted-map) (concat (map vector Bar Bar-y-B)
                                                         (map vector Brr Brr-y-B)
                                                         (map vector Baa Baa-y-B)))
                      y-AB (into (avl/sorted-map) (concat (map vector Aar Aar-y-AB)
                                                          (map vector Bar Bar-y-AB)
                                                          (map vector Aaa Aaa-y-AB)
                                                          (map vector Baa Baa-y-AB)))]
                  (assoc t
                         :Aar Aar
                         :Aaa Aaa
                         :Arr Arr
                         :Baa Baa
                         :Bar Bar
                         :Brr Brr
                         :gp-mean gp-mean-AB
                         :switching switching
                         :y {:A y-A :B y-B :AB y-AB}
                         :G G
                         :start-t start-t
                         :end-t end-t)))
        t-new (map t-mod trials)]
    (assoc AB :gp-var gp-var :gp-time-scale gp-time-scale :trials t-new)
    ))

(defn random-state [initial-state]
  (let [A (randomize-single-trials (:A initial-state))
        B (randomize-single-trials (:B initial-state))
        AB (randomize-dual-trials A B (:AB initial-state))]
    {:A A :B B :AB AB}))

(defn generate-dual-trials [num-trials gp-var gp-time-scale switching-prob A B]
  (let [{F-A :F gp-mean-A :gp-mean max-intensity-A :max-intensity} A
        {F-B :F gp-mean-B :gp-mean max-intensity-B :max-intensity} B
        f-A (comp first F-A)
        m-A (comp (partial + gp-mean-A) f-A)
        intensity-A (comp (partial * max-intensity-A) link m-A)
        f-B (comp first F-B)
        m-B (comp (partial + gp-mean-B ) f-B)
        intensity-B (comp (partial * max-intensity-B) link m-B)]
    {:gp-var gp-var
     :concentration 1
     :switch-p switching-prob
     :gp-time-locations [0.1 0.2 0.3 0.4 0.5]
     :gp-time-probabilities [0.2 0.2 0.2 0.2 0.2]
     :mean-prior-var 1
     :trials (take
              num-trials
              (repeatedly
               (fn [] (let [
                            switching (sample-bernoulli switching-prob)
                            ;switching (sample-binomial 1 :size 1 :prob switching-prob)
                            gp-time-scale (d/sample (d/uniform 0.1 0.11))
                            G (if (= switching 1) (sample-gp {} gp-var gp-time-scale) zero-3fn)
                            g (comp first G)
                            gp-mean-AB (d/sample (d/mixture [(d/normal -2 1) (d/normal 2 1)] [0.5 0.5]))
                            alpha (comp link (partial + gp-mean-AB) g)
                            [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr] (intensity-functions
                                                                   alpha
                                                                   intensity-A
                                                                   intensity-B
                                                                   max-intensity-A
                                                                   max-intensity-B)
                            start-t 0
                            end-t 1
                            Aaa (sample-ppp i-Aaa max-intensity-A start-t end-t)
                            Aar (sample-ppp i-Aar max-intensity-A start-t end-t)
                            Arr (sample-ppp i-Arr max-intensity-A start-t end-t)
                            Baa (sample-ppp i-Baa max-intensity-B start-t end-t)
                            Bar (sample-ppp i-Bar max-intensity-B start-t end-t)
                            Brr (sample-ppp i-Brr max-intensity-B start-t end-t)
                            obs-times (sort (concat Aaa Baa))
                            m-AB (comp (partial + gp-mean-AB) g)
                            ;Aar-y-A (map (fn [t] (sample-right-normal (m-A t) 1 )) Aar)
                            ;Aar-y-AB (map (fn [t] (sample-left-normal (m-AB t) 1 )) Aar)
                            ;Arr-y-A (map (fn [t] (sample-left-normal (m-A t) 1 )) Arr)
                            ;Bar-y-B (map (fn [t] (sample-right-normal (m-B t) 1 )) Bar)
                            ;Bar-y-AB (map (fn [t] (sample-right-normal (m-AB t) 1 )) Bar)
                            ;Brr-y-B (map (fn [t] (sample-left-normal (m-B t) 1 )) Brr)
                            ;Aaa-y-A (map (fn [t] (sample-right-normal (m-A t) 1 )) Aaa)
                            ;Aaa-y-AB (map (fn [t] (sample-right-normal (m-AB t) 1 )) Aaa)
                            ;Baa-y-B (map (fn [t] (sample-right-normal (m-B t) 1 )) Baa)
                            ;Baa-y-AB (map (fn [t] (sample-left-normal (m-AB t) 1 )) Baa)
                            Aar-y-A (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Aar)
                            Aar-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Aar)
                            Arr-y-A (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Arr)
                            Bar-y-B (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Bar)
                            Bar-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Bar)
                            Brr-y-B (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Brr)
                            Aaa-y-A (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Aaa)
                            Aaa-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Aaa)
                            Baa-y-B (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Baa)
                            Baa-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Baa)
                            y-A (into (avl/sorted-map) (concat (map vector Aar Aar-y-A)
                                                               (map vector Arr Arr-y-A)
                                                               (map vector Aaa Aaa-y-A)))
                            y-B (into (avl/sorted-map) (concat (map vector Bar Bar-y-B)
                                                               (map vector Brr Brr-y-B)
                                                               (map vector Baa Baa-y-B)))
                            y-AB (into (avl/sorted-map) (concat (map vector Aar Aar-y-AB)
                                                                (map vector Bar Bar-y-AB)
                                                                (map vector Aaa Aaa-y-AB)
                                                                (map vector Baa Baa-y-AB)))
                            ]
                        {:obs-times obs-times
                         :Aar Aar
                         :Aaa Aaa
                         :Arr Arr
                         :Baa Baa
                         :Bar Bar
                         :Brr Brr
                         :gp-mean gp-mean-AB
                         :dp-param gp-mean-AB
                         :dp-label 1
                         :gp-time-scale gp-time-scale
                         :switching switching
                         :y {:A y-A :B y-B :AB y-AB}
                         :G G
                         :start-t start-t
                         :end-t end-t}))))}))

(defn ulogpdf-gamma [a b x]
 (if (neg? x)
     Double/NEGATIVE_INFINITY
     (- (* (dec a) (log x)) (* b x))))
(defn logprior-gp-time-scale [x] (if (< x 0.09) (Double/NEGATIVE_INFINITY) (ulogpdf-gamma 3 20 x)))
(defn logprior-gp-var [x] (ulogpdf-gamma 10 2 x))

(defn update-dual-ys [m-A m-B trial]
  (let [{Aar :Aar
         Arr :Arr
         Bar :Bar
         Brr :Brr
         Baa :Baa
         Aaa :Aaa
         gp-mean-AB :gp-mean
         G :G} trial
        g (comp first G)
        m-AB (comp (partial + gp-mean-AB) g)
        Aar-y-A (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Aar)
        Aar-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Aar)
        Arr-y-A (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Arr)
        Bar-y-B (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Bar)
        Bar-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Bar)
        Brr-y-B (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Brr)
        Aaa-y-A (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Aaa)
        Aaa-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Aaa)
        Baa-y-B (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ 0.5 w) (/ 1 w)])) Baa)
        Baa-y-AB (map (fn [t] (let [w (d/sample (d/polya-gamma (m-AB t)))] [(/ -0.5 w) (/ 1 w)])) Baa)
        y-A (into (avl/sorted-map) (concat (map vector Aar Aar-y-A)
                                           (map vector Arr Arr-y-A)
                                           (map vector Aaa Aaa-y-A)))
        y-B (into (avl/sorted-map) (concat (map vector Bar Bar-y-B)
                                           (map vector Brr Brr-y-B)
                                           (map vector Baa Baa-y-B)))
        y-AB (into (avl/sorted-map) (concat (map vector Aar Aar-y-AB)
                                            (map vector Bar Bar-y-AB)
                                            (map vector Aaa Aaa-y-AB)
                                            (map vector Baa Baa-y-AB)))]
    (assoc trial :y {:A y-A :B y-B :AB y-AB})))

(defn update-dual-ys [m-A m-B trial]
  (let [{Aar :Aar
         Arr :Arr
         Bar :Bar
         Brr :Brr
         Baa :Baa
         Aaa :Aaa
         gp-mean-AB :gp-mean
         G :G} trial
        g (comp first G)
        m-AB (comp (partial + gp-mean-AB) g)
        Aar-y-A (map (fn [t] [(sample-right-normal (m-A t) 1 ) 1.0]) Aar)
        Aar-y-AB (map (fn [t] [(sample-left-normal (m-AB t) 1 ) 1.0]) Aar)
        Arr-y-A (map (fn [t] [(sample-left-normal (m-A t) 1 ) 1.0]) Arr)
        Bar-y-B (map (fn [t] [(sample-right-normal (m-B t) 1 ) 1.0]) Bar)
        Bar-y-AB (map (fn [t] [(sample-right-normal (m-AB t) 1 ) 1.0]) Bar)
        Brr-y-B (map (fn [t] [(sample-left-normal (m-B t) 1 ) 1.0]) Brr)
        Aaa-y-A (map (fn [t] [(sample-right-normal (m-A t) 1 ) 1.0]) Aaa)
        Aaa-y-AB (map (fn [t] [(sample-right-normal (m-AB t) 1 ) 1.0]) Aaa)
        Baa-y-B (map (fn [t] [(sample-right-normal (m-B t) 1 ) 1.0]) Baa)
        Baa-y-AB (map (fn [t] [(sample-left-normal (m-AB t) 1 ) 1.0]) Baa)
        y-A (into (avl/sorted-map) (concat (map vector Aar Aar-y-A)
                                           (map vector Arr Arr-y-A)
                                           (map vector Aaa Aaa-y-A)))
        y-B (into (avl/sorted-map) (concat (map vector Bar Bar-y-B)
                                           (map vector Brr Brr-y-B)
                                           (map vector Baa Baa-y-B)))
        y-AB (into (avl/sorted-map) (concat (map vector Aar Aar-y-AB)
                                            (map vector Bar Bar-y-AB)
                                            (map vector Aaa Aaa-y-AB)
                                            (map vector Baa Baa-y-AB)))]
    (assoc trial :y {:A y-A :B y-B :AB y-AB})))


(defn update-dual-times [intensity-A intensity-B max-intensity-A max-intensity-B trial]
  (let [{start-t :start-t
         end-t :end-t
         G :G
         gp-mean :gp-mean
         obs-times :obs-times} trial
        alpha (comp link (partial + gp-mean) first G)
        [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr] (intensity-functions
                                               alpha
                                               intensity-A
                                               intensity-B
                                               max-intensity-A
                                               max-intensity-B)
        Aar (sample-ppp i-Aar max-intensity-A start-t end-t)
        Arr (sample-ppp i-Arr max-intensity-A start-t end-t)
        Bar (sample-ppp i-Bar max-intensity-B start-t end-t)
        Brr (sample-ppp i-Brr max-intensity-B start-t end-t)
        AB-prob (fn [x] (/ (i-Aaa x) (+ (i-Aaa x) (i-Baa x))))
        AB-pred (fn [t] (= 1 (sample-bernoulli (AB-prob t))))
        ;AB-pred (fn [t] (= 1 (sample-binomial 1 :size 1 :prob (AB-prob t))))
        [Aaa Baa] (group-bool AB-pred obs-times)]
    (assoc trial :Aar Aar :Arr Arr :Bar Bar :Brr Brr :Baa Baa :Aaa Aaa)))

(defn update-dual-trials [state]
  (let [{F-A :F gp-mean-A :gp-mean max-intensity-A :max-intensity} (:A state)
        {F-B :F gp-mean-B :gp-mean max-intensity-B :max-intensity} (:B state)
        f-A (comp first F-A)
        m-A (comp (partial + gp-mean-A ) f-A)
        intensity-A (comp (partial * max-intensity-A) link m-A)
        f-B (comp first F-B)
        m-B (comp (partial + gp-mean-B ) f-B)
        intensity-B (comp (partial * max-intensity-B) link m-B)
        {trials :trials} (:AB state)
        update-times (partial update-dual-times intensity-A intensity-B max-intensity-A max-intensity-B)
        update-ys (partial update-dual-ys m-A m-B)
        update-trial (comp update-ys update-times)
        ]
    (assoc-in state [:AB :trials] (map update-trial trials))))

(defn dual-log-likelihood
  [gp-var gp-time-scale trial]
  (let [
        y (get-in trial [:y :AB])
        times (keys y)
        obs (map first (vals y))
        obs-var (map second (vals y))
        {gp-mean :gp-mean} trial
        centered-obs (map (fn [x] (- x gp-mean)) obs)
        ]
    (log-likelihood times centered-obs obs-var gp-var gp-time-scale)))


(defn log-sum-exp [coll]
  (let [a (reduce max coll)
        expshifted (map (fn [x] (exp (- x a))) coll)]
    (+ a (log (reduce + expshifted)))))

(defn switch-log-likelihood
  [switching gp-var gp-time-locations gp-time-probabilities trial]
  (let [{gp-time-scale :gp-time-scale} trial
        y (get-in trial [:y :AB])
        times (keys y)
        obs (map first (vals y))
        obs-var (map second (vals y))
        {gp-mean :gp-mean} trial
        ]
    (if (= switching 1)
      (let [centered-obs (map (fn [x] (- x gp-mean)) obs)
            loglik (map (fn [t] (log-likelihood times centered-obs obs-var gp-var t)) gp-time-locations)
            ;GP TIME LOCATIONS NOT ON LOG SCALE
            logprob (map log gp-time-probabilities)
            loglik+logprob (map + loglik logprob)]
        (log-sum-exp loglik+logprob))
      (reduce + 0 (map logpdf-normal obs (repeat gp-mean) obs-var)))))

(comment
  (defn switch-log-likelihood
    [switching gp-var trial]
    (let [{gp-time-scale :gp-time-scale} trial
          y (get-in trial [:y :AB])
          times (keys y)
          obs (map first (vals y))
          obs-var (map second (vals y))
          {gp-mean :gp-mean} trial
          ]
      (if (= switching 1)
        (log-likelihood times (map (fn [x] (- x gp-mean)) obs) obs-var gp-var gp-time-scale)
        (reduce + 0 (map logpdf-normal obs (repeat gp-mean) obs-var)))))

(defn switch-log-likelihood
  [switching gp-var trial]
  (let [{gp-time-scale :gp-time-scale} trial
        y (get-in trial [:y :AB])
        times (keys y)
        obs (map first (vals y))
        obs-var (map second (vals y))
        {gp-mean :gp-mean} trial
        ]
    (if (= switching 1)
      (log-likelihood times (map (fn [x] (- x gp-mean)) obs) obs-var gp-var gp-time-scale)
      (log-likelihood times (map (fn [x] (- x gp-mean)) obs) obs-var (* 0.2 gp-var) gp-time-scale))))

  (defn switch-log-likelihood
    [switching gp-var trial]
    (let [{gp-time-scale :gp-time-scale} trial
          y (get-in trial [:y :AB])
          times (keys y)
          obs (map first (vals y))
          obs-var (map second (vals y))
          ]
      (if (= switching 1)
        (log (d/quantile-integrate (comp exp (fn [mu] (log-likelihood times (map (fn [x] (- x mu)) obs) obs-var gp-var gp-time-scale))) (d/normal 0 s2-P) 20))
        (log (d/quantile-integrate (comp exp (fn [mu] (reduce + 0 (map logpdf-normal obs (repeat mu) obs-var)))) (d/normal 0 s2-P) 20))))))

(defn switch-probability [gp-var gp-time-locations gp-time-probabilities switch-p trial]
  (let [{gp-time-scale :gp-time-scale} trial
        on (+ (log switch-p) (switch-log-likelihood 1 gp-var gp-time-locations gp-time-probabilities trial))
        off (+ (log (- 1 switch-p)) (switch-log-likelihood 0 gp-var gp-time-locations gp-time-probabilities trial))
        odds (exp (- on off))]
    (/ odds (+ 1 odds))))

(defn dual-switching [gp-var gp-time-locations gp-time-probabilities switch-p trial]
  (let [{gp-time-scale :gp-time-scale} trial
        switching-probability (switch-probability gp-var gp-time-locations gp-time-probabilities switch-p trial)
        new-switching (if (< (rand) switching-probability) 1 0)]
    (assoc trial :switching new-switching)))

(defn update-dual-switch-probability [state]
  (let [trials (get-in state [:AB :trials])
        n (count trials)
        switchers (reduce + (map :switching trials))
        newprob (sample-beta 1 :alpha (+ 1 switchers) :beta (+ 1 (- n switchers)))]
    (assoc-in state [:AB :switch-p] newprob)))

(defn update-dual-switching [state]
  (let [{gp-var :gp-var
         switch-p :switch-p
         trials :trials
         gp-time-locations :gp-time-locations
         gp-time-probabilities :gp-time-probabilities} (:AB state)]
    (assoc-in state [:AB :trials] (map (partial dual-switching gp-var gp-time-locations gp-time-probabilities switch-p) trials))))

(defn gp-time-scale-posterior [loglikelihood locations probabilities]
  (let [
        logprobs (log probabilities)
        logliks (map loglikelihood locations)
        logpostprobs (map + logliks logprobs)
        maxprob (reduce max logpostprobs)
        offset-probs (map #(- % maxprob) logpostprobs)
        good-probs (exp offset-probs)
        Z (reduce + good-probs)
        newprobs (map #(/ % Z) good-probs)]
    (d/discrete-real locations good-probs)))


(defn dual-gp-time-scale [gp-var trial locations probabilities]
  (let [log-likelihood (if (= 1 (:switching trial)) (fn [t] (dual-log-likelihood gp-var t trial)) (fn [t] 0))]
    (assoc trial
           :gp-time-scale
           (d/sample (gp-time-scale-posterior log-likelihood locations probabilities)))))

(defn update-dual-gp-time-scale [{{gp-var :gp-var
                                   trials :trials
                                   locations :gp-time-locations
                                   probabilities :gp-time-probabilities
                                   } :AB :as state}]
  (let [new-trials (map #(dual-gp-time-scale gp-var % locations probabilities) trials)]
    (assoc-in state [:AB :trials] new-trials)))

(defn update-dual-gp-var [state]
  (let [{gp-var :gp-var
         trials :trials} (:AB state)
        active-trials (filter #(= 1 (:switching %)) trials)
        time-scales (map :gp-time-scale active-trials)
        conditional (fn [x]
                      (+ (logprior-gp-var x)
                         (reduce + 0
                                 (map (fn [y z] (dual-log-likelihood x y z)) time-scales active-trials))))]
    (assoc-in state [:AB :gp-var] (sample-slice conditional 1 gp-var))))


(defn dual-gp-mean [gp-var mean-prior-var trial]
  (if (= (:switching trial) 1)
    (let [{y :AB} (:y trial)
          {gp-time-scale :gp-time-scale} trial
          times (keys y)
          obs (map first (vals y))
          obs-var (map second (vals y))
          mus2 (mean-conditional times obs obs-var gp-var gp-time-scale)
          mu-L (first mus2)
          s2-L (second mus2)
          mu-P (:dp-param trial)
          posterior-var (/ 1  (+ (/ 1 s2-L) (/ 1 mean-prior-var)))
          posterior-mean (* posterior-var (+ (/ mu-P mean-prior-var) (/ mu-L s2-L)))
          new-mean (sample-normal 1 :mean posterior-mean :sd (sqrt posterior-var))]
      (assoc trial :gp-mean new-mean))
    (let [{y :AB} (:y trial)
          obs (map first (vals y))
          obs-var (map second (vals y))
          s2-L (/ 1 (reduce + 0 (div 1 obs-var)))
          mu-L (* s2-L (reduce + 0 (mul obs (div 1 obs-var))))
          mu-P (:dp-param trial)
          posterior-var (/ 1  (+ (/ 1 s2-L) (/ 1 mean-prior-var)))
          posterior-mean (* posterior-var (+ (/ mu-P mean-prior-var) (/ mu-L s2-L)))
          new-mean (sample-normal 1 :mean posterior-mean :sd (sqrt posterior-var)) ]
      (assoc trial :gp-mean new-mean))
    ))

(defn update-dual-gp-mean [state]
  (let [{gp-var :gp-var
         mean-prior-var :mean-prior-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)]
    (assoc-in state [:AB :trials] (map (partial dual-gp-mean gp-var mean-prior-var) trials))))

(defn dual-G [gp-var trial]
  (let [{switching :switching
         gp-time-scale :gp-time-scale} trial]
    (if (= switching 0)
      (assoc trial :G zero-3fn)
      (let [{y :AB} (:y trial)
            {gp-mean :gp-mean} trial
            times (keys y)
            obs (map first (vals y))
            obs-var (map second (vals y))
            centered-obs (map (fn [x] (- x gp-mean)) obs)
            pivots (FFBS times centered-obs obs-var gp-var gp-time-scale)]
        (assoc trial :G (sample-gp pivots gp-var gp-time-scale))))))

(defn update-dual-G [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)]
    (assoc-in state [:AB :trials] (map (partial dual-G gp-var) trials))))

(defn update-single-ys [gp-mean f trial]
  (let [{obs-times :obs-times
         aug-times :aug-times} trial
        obs-y (map (fn [t] [(sample-right-normal (+ gp-mean (f t)) 1 ) 1.0]) obs-times)
        aug-y (map (fn [t] [(sample-left-normal (+ gp-mean (f t)) 1 ) 1.0]) aug-times)
        ;obs-y (map (fn [t] (let [w (d/sample (d/polya-gamma (+ gp-mean (f t))))] [(/ 0.5 w) (/ 1 w)])) obs-times)
        ;aug-y (map (fn [t] (let [w (d/sample (d/polya-gamma (+ gp-mean (f t))))] [(/ -0.5 w) (/ 1 w)])) aug-times)
        obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
        aug-map (into (avl/sorted-map) (map vector aug-times aug-y))]
    (assoc trial :y (merge obs-map aug-map))))

(defn update-single-times [thin-intensity max-intensity trial]
  (let [{start-t :start-t
         end-t :end-t} trial]
    (assoc trial :aug-times (sample-ppp thin-intensity max-intensity start-t end-t))))

(defn update-single-trials [typ state]
  (let [{F :F
         max-intensity :max-intensity
         gp-mean :gp-mean
         trials :trials} (typ state)
        f (comp first F)
        intensity (comp (partial * max-intensity) link (partial + gp-mean) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))
        update-ys (partial update-single-ys gp-mean f)
        update-single-times (partial update-single-times thin-intensity max-intensity)
        update-trial (comp update-ys update-single-times)
        ]
    (assoc-in state [typ :trials] (map update-trial trials))))

(defn update-single-gp-mean [times obs obs-var typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale} (typ state)
        mus2 (mean-conditional times obs obs-var gp-var gp-time-scale)
        new-mean (sample-normal 1 :mean (first mus2) :sd (sqrt (second mus2)))]
    (assoc-in state [typ :gp-mean] new-mean)))

(defn update-single-gp-time-scale-approx [typ state]
  (let [{gp-var :gp-var
         F :F} (typ state)
        logconditional (fn [x] (logpdf-prior-gp-approx F gp-var x))
        logliks (map logconditional [0.1 0.2 0.3 0.4 0.5])]
    (assoc-in state [typ :gp-time-scale] (d/sample (d/discrete-real [0.1 0.2 0.3 0.4 0.5] logliks :log? true)))))

(defn update-single-gp-time-scale [times obs obs-var typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         gp-mean :gp-mean} (typ state)
        centered-obs (map (fn [x] (- x gp-mean)) obs)
        conditional (fn [x] (+ (logprior-gp-time-scale x)
                               (log-likelihood times centered-obs obs-var gp-var x)))]
    (assoc-in state [typ :gp-time-scale] (sample-slice conditional 1 gp-time-scale))))

(defn update-single-gp-var-approx [typ state]
  (let [{F :F
         gp-time-scale :gp-time-scale} (typ state)]
    (assoc-in state [typ :gp-var] (d/sample (variance-full-conditional-approx F gp-time-scale)))))

(defn update-single-gp-var [times obs obs-var typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         gp-mean :gp-mean} (typ state)
        centered-obs (map (fn [x] (- x gp-mean)) obs)
        conditional (fn [x] (+ (logprior-gp-var x)
                               (log-likelihood times centered-obs obs-var x gp-time-scale)))]
    (assoc-in state [typ :gp-var] (sample-slice conditional 1 gp-var))))

(defn update-single-F [times obs obs-var typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         gp-mean :gp-mean} (typ state)
        centered-obs (map (fn [x] (- x gp-mean)) obs)
        pivots (FFBS times centered-obs obs-var gp-var gp-time-scale)]
    (assoc-in state [typ :F] (sample-gp pivots gp-var gp-time-scale))))

(defn update-single-intensity-approx [typ state]
  (let [single-ys (map :y (get-in state [typ :trials]))
        dual-ys (map (fn [x] (get-in x [:y typ])) (get-in state [:AB :trials]))
        y (apply merge (concat single-ys dual-ys))
        times (keys y)
        obs (map first (vals y))
        obs-var (map second (vals y))
        update-gp-mean (partial update-single-gp-mean times obs obs-var typ)
        update-gp-time-scale (partial update-single-gp-time-scale-approx typ)
        update-gp-var (partial update-single-gp-var-approx typ)
        update-F (partial update-single-F times obs obs-var typ)
        update (comp update-gp-var
                     update-gp-time-scale
                     update-F
                     update-gp-mean
                     )]
    (update state)))

(defn update-single-intensity [typ state]
  (let [single-ys (map :y (get-in state [typ :trials]))
        dual-ys (map (fn [x] (get-in x [:y typ])) (get-in state [:AB :trials]))
        y (apply merge (concat single-ys dual-ys))
        times (keys y)
        obs (map first (vals y))
        obs-var (map second (vals y))
        update-gp-mean (partial update-single-gp-mean times obs obs-var typ)
        update-gp-time-scale (partial update-single-gp-time-scale times obs obs-var typ)
        update-gp-var (partial update-single-gp-var times obs obs-var typ)
        update-F (partial update-single-F times obs obs-var typ)
        update (comp update-F
                     update-gp-var
                     update-gp-time-scale
                     update-gp-mean
                     )]
    (update state)))

(def dirichlet-alpha [0.2 0.2 0.2 0.2 0.2])

(defn update-dp-labels [state]
  (let [{concentration :concentration
         mean-prior-var :mean-prior-var
         trials :trials} (:AB state)
        labels (map :dp-label trials)
        observations (map :gp-mean trials)
        newlabels (d/dp-update-labels observations (d/normal :mu mean-prior-var) base-measure concentration (into [] labels))
        newparams (d/dp-update-params observations (d/normal :mu mean-prior-var) base-measure concentration newlabels)]
    (assoc-in state [:AB :trials] (map (fn [x y z] (assoc x :dp-label y :dp-param z)) trials newlabels newparams))
    ))

(defn update-probabilities [state]
  (let [{probabilities :gp-time-probabilities
         locations :gp-time-locations
         trials :trials} (:AB state)
        n (for [t locations] (reduce + (map (fn [x] (if (= t (:gp-time-scale x)) 1 0)) trials)))
        new-alpha (add n dirichlet-alpha)]
    (assoc-in state [:AB :gp-time-probabilities] (d/sample (d/dirichlet new-alpha)))))

(defn update-mean-prior-var [state]
  (let [{trials :trials} (:AB state)
        means (map :gp-mean trials)
        dp-params (map :dp-param trials)
        reducer (fn [p [m n]] (d/posterior [m] (d/normal n :s2) p))
        post (reduce reducer mean-prior-var-prior (map vector means dp-params))]
    (assoc-in state [:AB :mean-prior-var] (d/sample post))))

(def transition (comp
                 (partial update-single-intensity :B)
                 (partial update-single-trials :B)
                 (partial update-single-intensity :A)
                 (partial update-single-trials :A)
                 update-dual-switch-probability
                 update-mean-prior-var
                 update-dp-labels
                 update-dual-G
                 update-probabilities
                 update-dual-gp-var
                 update-dual-gp-mean
                 update-dual-gp-time-scale
                 update-dual-switching
                 update-dual-trials
                 )
  )

(defn strip-data [coll]
  (prewalk #(if (map? %) (dissoc % :y :Baa :Bar :Brr :Aaa :Aar :Arr :obs-times :aug-times) %) coll))

(defn strip-intensity [coll]
  (prewalk #(if (clojure.test/function? %) (map (comp first %) (range 0 1 0.001)) %) coll))

(defn functionify [coll]
  (let [gridded (into (avl/sorted-map) (map vector (range 0 1 (/ 1 (count coll))) coll))]
    (fn [x] [(second (avl/nearest gridded <=  x)) 0 0])))

(defn functionify-intensity [coll]
  (prewalk #(if (sequential? %) (if (= (count %) 1000) (functionify %)  %) %) coll))



(def strip (comp strip-intensity strip-data))
(def strip strip-data)

(defn take-thin [take thin transition initial-state]
  (loop [iter 0
         acc '()
         acclen 0
         state initial-state]
    (if (= take acclen)
      acc
      (if (zero? (mod iter thin))
        (do
          (println iter)
          (recur (inc iter) (conj acc (strip state)) (inc acclen) (transition state)))
        (recur (inc iter) acc acclen (transition state))))))


(def Atrials (generate-single-trials 20 2 0.1 0.1 400))
(def Btrials (generate-single-trials 20 2 0.1 0.1 50))
(def ABtrials (generate-dual-trials 20 5 0.1 0.5 Atrials Btrials))
(def initial-state {:A Atrials :B Btrials :AB ABtrials})
;(def mcmc (iterate transition initial-state))
;(def mcmc (iterate transition (random-state initial-state)))
;(:gp-time-scale (:AB (assoc-in initial-state [:AB :gp-time-scale] 0.2)))

;(def iterates (take-thin 1000 10 transition (assoc-in initial-state [:AB :gp-time-scale] 0.3)))
(def run-time (with-out-str (time (def iterates (take-thin 1000 10 transition initial-state)))))
(ic/view (histogram (sample-normal 1000)))
(def iterates2 (take-thin 1000 10 transition initial-state))
(def foo (doall (update-dual-G (update-probabilities (update-dual-gp-time-scale (update-dual-gp-var (update-dual-gp-mean (update-dual-switching (update-dual-trials initial-state)))))))))
(last iterates)
(def iterates (take 2 (take-nth 2000 mcmc)))
(def iterates (take-thin 200 10 mcmc))
(last iterates)
(def iterates2 1)
(def mcmc 1)
(System/gc)
(:gp-time-scale (second (:trials ABtrials)))
(:trials ABtrials)
(first (:trials  ABtrials))
(:A (:y (update-dual-ys-pg (fn [x] 0) (fn [x] 0) (first (:trials ABtrials)))))
(spit "initial-state" (json/write-str (strip-intensity initial-state)))
(spit "probit-surya-data" (json/write-str iterates2))
(def surya-data (json/read-str (slurp "surya-data")))
(defn convert-to-seconds [x] (map #(/ % 1000) x))
(def Atimes (map convert-to-seconds (get-in surya-data ["spiketimes" "A"])))
(def Btimes (map convert-to-seconds (get-in surya-data ["spiketimes" "B"])))
(def ABtimes (map convert-to-seconds (get-in surya-data ["spiketimes" "AB"])))
(def Aemptytrials (generate-single-trials (count Atimes) 2 0.1 0.1 400))
(def Bemptytrials (generate-single-trials (count Btimes) 2 0.1 0.1 100))
(def Atrials (assoc Aemptytrials :trials (map (fn [[trial times]] (assoc trial :obs-times times)) (map vector (:trials Aemptytrials) Atimes))))
(def Btrials (assoc Bemptytrials :trials (map (fn [[trial times]] (assoc trial :obs-times times)) (map vector (:trials Bemptytrials) Btimes))))
(def ABemptytrials (generate-dual-trials (count ABtimes) 5 0.1 0.5 Atrials Btrials))
(def ABtrials (assoc ABemptytrials :trials (map (fn [[trial times]] (assoc trial :obs-times times)) (map vector (:trials ABemptytrials) ABtimes))))
(def initial-state {:A Atrials :B Btrials :AB ABtrials})

(def time-transition (comp
                      (time (partial update-single-intensity :B))
                      (time (partial update-single-trials :B))
                      (time (partial update-single-intensity :A))
                      (time (partial update-single-trials :A))
                      (time update-dual-switch-probability)
                      (time update-mean-prior-var)
                      (time update-dp-labels)
                      (time update-dual-G)
                      (time update-probabilities)
                      (time update-dual-gp-var)
                      (time update-dual-gp-mean)
                      (time update-dual-gp-time-scale)
                      (time update-dual-switching)
                      (time update-dual-trials)
                 )
  )
(def test-state (doall (last (take 3 (iterate transition initial-state)))))
(def foo (time (transition test-state)))
(def foo (time (doall (update-dual-trials test-state))))
(def foo (time (doall (update-dual-switching test-state))))
(def foo (time (doall (update-dual-gp-time-scale test-state))))
(def foo (time (doall (update-dual-gp-mean test-state))))
(def foo (time (doall (update-dual-gp-var test-state))))
(def foo (time (doall (update-probabilities test-state))))
(def foo (time (doall (update-dual-G test-state))))
(def foo (time (doall (update-dp-labels test-state))))
(def foo (time (doall (update-mean-prior-var test-state))))
(def foo (time (doall (update-dual-switch-probability test-state))))
(def foo (time (doall (update-single-trials :A test-state))))
(def foo (time (doall (update-single-intensity :A test-state))))
(def foo (time (doall (update-single-trials :B test-state))))
(def foo (time (doall (update-single-intensity :B test-state))))
(def single-ys (map :y (get-in test-state [:A :trials])))
(def dual-ys (map (fn [x] (get-in x [:y :A])) (get-in test-state [:AB :trials])))
(def y (take 100000 (apply merge (concat single-ys dual-ys))))
(def times (keys y))
(def obs (map first (vals y)))
(def obs-var (map second (vals y)))
(def update-gp-mean (partial update-single-gp-mean times obs obs-var :A))
(def update-gp-time-scale (partial update-single-gp-time-scale times obs obs-var :A))
(def update-gp-var (partial update-single-gp-var times obs obs-var :A))
(def update-gp-time-scale-approx (partial update-single-gp-time-scale-approx :A))
(def update-gp-var-approx (partial update-single-gp-var-approx :A))
(def update-F (partial update-single-F times obs obs-var :A))
(def foo (time (update-gp-mean test-state)))
(def foo (time (update-gp-time-scale test-state)))
(def foo (time (update-gp-var test-state)))
(def foo (time (update-F test-state)))
(def foo (time (update-gp-time-scale-approx test-state)))
(def foo (time (update-gp-var-approx test-state)))
(count @(:data (meta (:F (:A test-state)))))
(count times)
()

(def new-state (first (take 3 (iterate transition initial-state))))

(keys (first iterates))
(def iterates (json/read-str (slurp "probit-surya-data") :key-fn keyword))
(def iterates2 (map functionify-intensity  iterates))
(count iterates)
(= 1000 (count (get-in (functionify-intensity (first iterates2)) [:A :F])))
(seq? (get-in (functionify-intensity (first iterates2)) [:A :F]))
(get-in (first iterates2) [:A :F])
(+ 1 2)
(get-in (first (map functionify-intensity iterates)) [:A :F])
(reduce (fn [x y] (and x y)) (map #(= (type %) clojure.lang.PersistentVector) (map second (:A (:y (update-dual-ys-link (fn [x] 2) (fn [x] 3) (first (:trials ABtrials)))))) ))
(map second (vals (:y (first (:trials (:B (nth iterates 10)))))))
(map second (vals (:y (first (:trials (:A (nth iterates 10)))))))
(map second (vals (:A (:y (first (:trials (:AB (nth iterates 10))))))))
(map second (vals (:B (:y (first (:trials (:AB (nth iterates 10))))))))
(map second (vals (:AB (:y (first (:trials (:AB (nth iterates 10))))))))


(first (:F (:A (last iterates2))))



(ic/save (histogram (map (fn [x] (get-in x [:AB :mean-prior-var])) iterates) :nbins 100) "figures/mean-prior-var.png")
(ic/save (histogram (map (fn [x] (get-in x [:AB :switch-p])) iterates) :nbins 100) "figures/cell-switch-probabiity.png")
(defn plot-dp-density [iterates]
  (let [concentrations (map (fn [x] (get-in x [:AB :concentration])) iterates)
        mean-prior-vars (map (fn [x] (get-in x [:AB :mean-prior-var])) iterates)
        trials (map (fn [x] (get-in x [:AB :trials])) iterates)
        params (map (fn [x] (map :dp-param x)) trials)
        pred-fn-draws (map (fn [params mean-prior-var concentration] (d/pdf (d/marginal (d/normal :mu mean-prior-var) (d/posterior-predictive params :G (d/dirichlet-process concentration base-measure))))) params mean-prior-vars concentrations )
        lower-function (fn [x] (quantile (map (fn [y] (y x)) pred-fn-draws) :probs 0.025))
        upper-function (fn [x] (quantile (map (fn [y] (y x)) pred-fn-draws) :probs 0.975))
        mean-function (fn [x] (mean (map (fn [y] (y x)) pred-fn-draws)))
        ]
    (doto
        (function-plot mean-function -5 5)
      (add-function lower-function -5 5)
      (add-function upper-function -5 5)
      (set-stroke :dash 5 :series 0) 
      )

    ))

(ic/view (plot-dp-density iterates))


;;;;;;;;;;;;;;;;;;GORILLA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn plot-mean-prior-var [iterates]
  (plot/histogram (map (fn [x] (get-in x [:AB :mean-prior-var])) iterates) :bins 100))

(defn plot-switch-probability [iterates]
  (plot/histogram (map (fn [x] (get-in x [:AB :switch-p])) iterates) :bins 100))

(defn plot-dp-density [iterates]
  (let [concentration (:concentration (:AB (first iterates)))
        mean-prior-vars (map (fn [x] (get-in x [:AB :mean-prior-var])) iterates)
        trials (map (fn [x] (get-in x [:AB :trials])) iterates)
        params (map (fn [x] (map :dp-param x)) trials)
        pred-fn-draws (map (fn [params mean-prior-var] (d/pdf (d/marginal (d/normal :mu mean-prior-var) (d/posterior-predictive params :G (d/dirichlet-process concentration base-measure))))) params mean-prior-vars)
        ]
    (apply plot/compose (map (fn [p] (plot/plot p [-5 5] :colour "skyblue" :opacity 0.1)) pred-fn-draws))
    ))


(defn plot-intensity
  ([typ iterates]
   (let[
        ;mean-function (fn [x] (mean (map (fn [y] (* (:max-intensity y) (link (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates))))
        lower-function (fn [x] (quantile (map (fn [y] (* (:max-intensity y) (link (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates)) :probs 0.025))
        upper-function (fn [x] (quantile (map (fn [y] (* (:max-intensity y) (link (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates)) :probs 0.975))
        ]
     (plot/compose (plot/plot upper-function [0 1] :plot-range [[0 1] [0 410]] :colour "skyblue")
                   (plot/plot lower-function [0 1] :plot-range [[0 1] [0 410]] :colour "skyblue"))
     ))
  ([typ iterates truth]
   (let [true-function (fn [x] (* (:max-intensity truth) (link (+ (:gp-mean truth) ((comp first (:F truth)) x)))))]
     (plot/compose
      (plot-intensity typ iterates)
      (plot/plot true-function [0 1] :plot-range [[0 1] [0 410]] :colour "#FA8072"))
     )
   ))

(defn plot-gp
  ([typ iterates]
   (let[
        ;mean-function (fn [x] (mean (map (fn [y] (* (:max-intensity y) (link (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates))))
        lower-function (fn [x] (quantile (map (fn [y] (first ((:F y) x))) (map typ iterates)) :probs 0.025))
        upper-function (fn [x] (quantile (map (fn [y] (first ((:F y) x))) (map typ iterates)) :probs 0.975))
        ]
     (plot/compose (plot/plot upper-function [0 1] :plot-range [[0 1] [-4 4]] :colour "skyblue")
                   (plot/plot lower-function [0 1] :colour "skyblue"))
     ))
  ([typ iterates truth]
   (let [true-function (fn [x] ((comp first (:F truth)) x))]
     (plot/compose
      (plot-gp typ iterates)
      (plot/plot true-function [0 1] :colour "#FA8072"))
     )
   ))


(defn plot-param [param typ iterates]
  (plot/histogram (map (fn [x] (param (typ x))) iterates) :normalise :probability-density :bins 30))

(defn plot-trace [param typ iterates]
  (plot/list-plot (map (fn [x] (param (typ x))) iterates) :joined true) )

(defn quartile-functions
  [function-stream]
  (let [lower-function (fn [x] (quantile (map (fn [y] (y x)) function-stream) :probs 0.025))
        upper-function (fn [x] (quantile (map (fn [y] (y x)) function-stream) :probs 0.975))]
    [lower-function upper-function]))

(defn extract-params [iterates]
  (let [max-intensity-A (map #(get-in % [:A :max-intensity]) iterates)
        max-intensity-B (map #(get-in % [:B :max-intensity]) iterates)
        f-A (map #(comp first (get-in % [:A :F])) iterates)
        gp-mean-A (map #(get-in % [:A :gp-mean]) iterates)
        f-B (map #(comp first (get-in % [:B :F])) iterates)
        gp-mean-B (map #(get-in % [:B :gp-mean]) iterates)
        intensity-A (map (fn [Lambda gp-mean f] (fn [t] (* Lambda (link (+ gp-mean (f t)))))) max-intensity-A gp-mean-A f-A)
        intensity-B (map (fn [Lambda gp-mean f] (fn [t] (* Lambda (link (+ gp-mean (f t)))))) max-intensity-B gp-mean-B f-B)]
    {:intensity-A intensity-A :intensity-B intensity-B}
    ))

(defn plot-gp-dual
  ([g-AB]
   (let [[lower-g upper-g] (quartile-functions g-AB)
         g-plot (plot/compose (plot/plot upper-g [0 1] :plot-range [[0 1] [-4 4]] :colour "skyblue")
                              (plot/plot lower-g [0 1] :colour "skyblue"))]
     g-plot))
  ([g-AB i truth]
   (let [g-true (comp first (:G (nth (get-in truth [:AB :trials]) i)))]
     (plot/compose
      (plot/plot g-true [0 1] :plot-range [[0 1] [-4 4]] :colour "#FA8072")
      (plot-gp-dual g-AB)))
   ))

(defn get-truth [i truth]
  (let [max-intensity-A (get-in truth [:A :max-intensity])
        max-intensity-B (get-in truth [:B :max-intensity])
        f-A (comp first (get-in truth [:A :F]))
        gp-mean-A (get-in truth [:A :gp-mean])
        f-B (comp first (get-in truth [:B :F]))
        gp-mean-B (get-in truth [:B :gp-mean])
        g-AB (comp first (:G (nth (get-in truth [:AB :trials]) i)))
        gp-mean-AB (:gp-mean (nth (get-in truth [:AB :trials]) i))
        intensity-A (fn [t] (* max-intensity-A (link (+ gp-mean-A (f-A t)))))
        intensity-B (fn [t] (* max-intensity-B (link (+ gp-mean-B (f-B t)))))
        alpha (fn [t] (link (+ gp-mean-AB (g-AB t))))
        intensity-AB (fn [t] (+ (* (alpha t) (intensity-A t))
                                (* (- 1 (alpha t)) (intensity-B t))))]
    {:intensity-AB intensity-AB
     :g-AB g-AB}
    ))

(defn plot-intensity-dual
  ([i iterates gp-mean-AB g-AB]
   (let [{intensity-A :intensity-A
          intensity-B :intensity-B} (extract-params iterates)
         alpha (map (fn [gp-mean g] (fn [t] (link (+ gp-mean (g t))))) gp-mean-AB g-AB)
         intensity-AB (map (fn [in-A in-B a]
                             (fn [t] (+ (* (a t) (in-A t))
                                        (* (- 1 (a t)) (in-B t))))) intensity-A intensity-B alpha)
         [lower-intensity upper-intensity] (quartile-functions intensity-AB)]
     (plot/compose (plot/plot upper-intensity [0 1] :plot-range [[0 1] [0 410]] :colour "skyblue")
                   (plot/plot lower-intensity [0 1] :colour "skyblue"))
     ))
  ([i iterates truth gp-mean-AB g-AB]
   (let [{intensity-AB :intensity-AB} (get-truth i truth)]
     (plot/compose
      (plot/plot intensity-AB [0 1] :plot-range [[0 1] [0 410]] :colour "#FA8072")
      (plot-intensity-dual i iterates gp-mean-AB g-AB)))
   ))


(defn plot-dual [i iterates truth]
  (let [g-AB (map #(comp first (:G (nth (get-in % [:AB :trials]) i))) iterates)
        true-gp-mean-AB (:gp-mean (nth (get-in truth [:AB :trials]) i))
        gp-mean-AB (map #(:gp-mean (nth (get-in % [:AB :trials]) i)) iterates)
        gp-time-scale (map #(:gp-time-scale (nth (get-in % [:AB :trials]) i)) iterates)
        switching (map #(:switching (nth (get-in % [:AB :trials]) i)) iterates)
        switch-prob (mean switching)
        switch-trace (plot/list-plot switching :joined true)
        gp-time-scale-trace (plot/histogram gp-time-scale :bins 30)
        switch-plot (plot/bar-chart ["Switching" "Non Switching"] [switch-prob (- 1 switch-prob)])
        mean-histogram (plot/histogram gp-mean-AB :bins 30)
        mean-traceplot (plot/list-plot gp-mean-AB :joined true)
        g-plot (plot-gp-dual g-AB i truth)
        intensity-plot (plot-intensity-dual i iterates truth gp-mean-AB g-AB)
        ]
    [true-gp-mean-AB mean-histogram mean-traceplot gp-time-scale-trace switch-trace switch-plot g-plot intensity-plot]))

(defn plot-duals [iterates truth]
  (map #(plot-dual % iterates truth) (range 0 (count (get-in (first iterates) [:AB :trials])))))

(defn plot-probs [iterates]
  (let [probs (map (fn [x] (:gp-time-probabilities (:AB x))) iterates)
        summedprobs (reduce add probs)
        normalprobs (div summedprobs (count iterates))]
    normalprobs))


 (get-in (first iterates) [:AB :trials])
(defn plot-gp-duals [iterates true-Gs]
  (let [Gstreams (for [i (range 0 (count (get-in (first iterates) [:AB :trials])))]
                   (map #(:G (nth (get-in % [:AB :trials]) i)) iterates))]
    (map plot-gp-dual Gstreams true-Gs)))

(use '[gorilla-repl latex table])
(use 'gorilla-renderable.core)

(defn summarize [iterates]
 (table-view 
  [[(latex-view "\\lambda_A(t)=\\Lambda \\Phi (\\mu_A+f_A(t))") 
    "Intensity Function of A-trials. Truth red" 
    (plot-intensity :A iterates Atrials)
    "NA"]
   [(latex-view "f_A(t)") 
    "GP Function of A-trials. Truth red" 
    (plot-gp :A iterates Atrials)
    "NA"]
   [(latex-view "\\mu_A") 
    (str "Defines mean Intensity of A-trials. Truth " (:gp-mean Atrials)) 
    (plot-param :gp-mean :A iterates)
    (plot-trace :gp-mean :A iterates)]
   [(latex-view "\\tau_A") 
    (str "Characteristic time-scale of Gaussian Process f. Truth " (:gp-time-scale Atrials)) 
    (plot-param :gp-time-scale :A iterates)
    (plot-trace :gp-time-scale :A iterates)]
   [(latex-view "\\sigma^2_A") 
    (str "Variance of Gaussian Process f. Truth " (:gp-var Atrials)) 
    (plot-param :gp-var :A iterates)
    (plot-trace :gp-var :A iterates)
    ]
   [(latex-view "\\lambda_B(t)=\\Lambda \\Phi (\\mu_B+f_B(t))") 
    "Intensity Function of B-trials" 
    (plot-intensity :B iterates Btrials)
    "NA"]
   [(latex-view "f_B(t)") 
    "GP Function of B-trials. Truth red" 
    (plot-gp :B iterates Btrials)
    "NA"]
   [(latex-view "\\mu_B") 
    (str "Defines mean Intensity of B-trials. Truth " (:gp-mean Btrials)) 
    (plot-param :gp-mean :B iterates)
    (plot-trace :gp-mean :B iterates)]
    [(latex-view "\\tau_B") 
    (str "Characteristic time-scale of Gaussian Process f_B. Truth " (:gp-time-scale Btrials)) 
    (plot-param :gp-time-scale :B iterates)
    (plot-trace :gp-time-scale :B iterates)]
   [(latex-view "\\sigma^2_B") 
    (str "Variance of Gaussian Process f_B. Truth " (:gp-var Btrials)) 
    (plot-param :gp-var :B iterates)
    (plot-trace :gp-var :B iterates)]
   [(latex-view "\\sigma^2_{AB}") 
    (str "Variance of all Gaussian Processes g_AB. Truth " (:gp-var ABtrials)) 
    (plot-param :gp-var :AB iterates)
    (plot-trace :gp-var :AB iterates)] ]
    :columns 
  ["Paramter" "Description" "Posterior Summary" "TracePlot"]) )

(def gp-var (:gp-var ABtrials))
(def gp-time-scale (:gp-time-scale ABtrials))
(def m-A (fn [x] (+ (:gp-mean Atrials) (first ((:F Atrials) x)))))
(def m-B (fn [x] (+ (:gp-mean Btrials) (first ((:F Btrials) x)))))
(def intensity-A (fn [x] (* 100 (link (m-A x)))))
(def intensity-B (fn [x] (* 100 (link (m-B x)))))
(require '[ssm4clj.gp :refer :all])
(def update-times (partial update-dual-times intensity-A intensity-B 100))
(def update-ys (partial update-dual-ys m-A m-B))
(def update-ys-pg (partial update-dual-ys-pg m-A m-B))
(def update-G (partial dual-G gp-var gp-time-scale))

(def trial (nth (:trials ABtrials 4)))
(plot/plot (comp first (:G trial)) [0 1])
(switch-probability gp-var gp-time-scale trial)
(def midtrial (assoc trial :switching 0))
(def newtrial (update-ys (update-times (update-G midtrial))))
(switch-probability gp-var gp-time-scale newtrial)

(defn rpg
  ([n b c]
   (take n (repeatedly #(rpg b c))))
  ([b c]
   (/ (reduce + (map (fn [k] (/  (sample-gamma 1 :shape b :scale 1) (+ (square (- k 0.5)) (/ (square c) (* 4 (square Math/PI)))))) (range 1 10000))) (* 2 (square Math/PI)))))

(rpg 100 4 3)
(sample-gamma :shape 10 :scale 0.3)


(zero-matrix 2 2)
(count t)
(fill (zero-matrix (count t) (count t)) 1)

(:gp-mean trial)
(d/mvnor)




(def trial (nth (:trials ABtrials) 3))
(plot/plot (comp first (:G trial)) [0 1])
(plot/list-plot (:AB (:y trial)))
(switch-probability gp-var gp-time-scale trial)
(def midtrial (assoc trial :switching 0))
(def newtrial (update-ys (update-times (update-G midtrial))))
(plot/list-plot (:AB (:y newtrial)))
(switch-probability gp-var gp-time-scale newtrial)


(def trial (update-ys-pg (update-times (nth (:trials ABtrials) 3))))
(plot/plot (comp first (:G trial)) [0 1])
(def t (keys (:AB (:y trial))))
(def y (map first (vals (:AB (:y trial)))))
(def omegas (map second (vals (:AB (:y trial)))))
(plot/list-plot (map vector t y ))
(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (diagonal-matrix (div 1 omegas)))
(def covmat (add K obsvars))
(def meanvector (fill (zero-vector (count t)) (:gp-mean trial)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector obsvars) y))
(def odds (exp (- loglik1 loglik0)))
(/ odds (+ 1 odds ))
;(switch-probability gp-var gp-time-scale trial)
(def midtrial (assoc trial :switching 0))
(def newtrial (update-ys-pg (update-times (update-G midtrial))))
(plot/plot (comp first (:G newtrial)) [0 1])
(def t (keys (:AB (:y newtrial))))
(def y (map first (vals (:AB (:y newtrial)))))
(def omegas (map second (vals (:AB (:y newtrial)))))
(plot/list-plot (map vector t y ))
(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (diagonal-matrix (div 1 omegas)))
(def covmat (add K obsvars))
(def meanvector (fill (zero-vector (count t)) (:gp-mean newtrial)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector obsvars) y))
(def odds (exp (- loglik1 loglik0)))
(/ odds (+ 1 odds ))
;(plot/list-plot (:AB (:y newtrial)))
;(switch-probability gp-var gp-time-scale newtrial)



(def polyahist (take 200 (repeatedly #(do (def trial (update-ys-pg (update-times (nth (:trials ABtrials) 3))))

(def t (keys (:AB (:y trial))))
(def y (map first (vals (:AB (:y trial)))))
(def omegas (map second (vals (:AB (:y trial)))))

(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (diagonal-matrix (div 1 omegas)))
(def covmat (add K obsvars))
(def meanvector (fill (zero-vector (count t)) (:gp-mean trial)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector obsvars) y))
(def odds (exp (- loglik1 loglik0)))

;(switch-probability gp-var gp-time-scale trial)
(def midtrial (assoc trial :switching 0))
(def newtrial (update-ys-pg (update-times (update-G midtrial))))

(def t (keys (:AB (:y newtrial))))
(def y (map first (vals (:AB (:y newtrial)))))
(def omegas (map second (vals (:AB (:y newtrial)))))

(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (diagonal-matrix (div 1 omegas)))
(def covmat (add K obsvars))
(def meanvector (fill (zero-vector (count t)) (:gp-mean newtrial)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector obsvars) y))
(def odds (exp (- loglik1 loglik0)))
(/ odds (+ 1 odds ))))))
;(plot/list-plot (:AB (:y newtrial)))
;(switch-probability gp-var gp-time-scale newtrial)

(def linkhist (take 200 (repeatedly #(do (def trial (nth (:trials ABtrials) 3))
                                           (def midtrial (assoc trial :switching 0))
                                           (def newtrial (update-ys (update-times (update-G midtrial))))
                                           (switch-probability gp-var gp-time-scale newtrial)))))


(def linkhist2 (take 200 (repeatedly #(do (def trial (update-ys (update-times (nth (:trials ABtrials) 3))))

(def t (keys (:AB (:y trial))))
(def y  (vals (:AB (:y trial))))

(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (identity-matrix (count t)))
(def covmat (add K obsvars))
(def meanvector (fill (zero-vector (count t)) (:gp-mean trial)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector obsvars) y))
(def odds (exp (- loglik1 loglik0)))

;(switch-probability gp-var gp-time-scale trial)
(def midtrial (assoc trial :switching 0))
(def newtrial (update-ys (update-times (update-G midtrial))))

(def t (keys (:AB (:y newtrial))))
(def y (vals (:AB (:y newtrial))))

(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (identity-matrix (count t)))
(def covmat (add K obsvars))
(def meanvector (fill (zero-vector (count t)) (:gp-mean newtrial)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector obsvars) y))
(def odds (exp (- loglik1 loglik0)))
(/ odds (+ 1 odds ))))))
;(plot/list-plot (:AB (:y newtrial)))
;(switch-probability gp-var gp-time-scale newtrial)


(def polyahist2 (take 200 (repeatedly #(do (def trial (update-ys-pg (update-times (nth (:trials ABtrials) 3))))

(def t (keys (:AB (:y trial))))
(def y (map first (vals (:AB (:y trial)))))
(def omegas (map second (vals (:AB (:y trial)))))

(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (diagonal-matrix (div 1 omegas)))
(def covmat (add K obsvars (fill (zero-matrix (count t) (count t)) (:gp-mean trial))))
(def meanvector (zero-vector (count t)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector (add obsvars (fill (zero-matrix (count t) (count t)) 1))) y))
(def odds (exp (- loglik1 loglik0)))

;(switch-probability gp-var gp-time-scale trial)
(def midtrial (assoc trial :switching 0))
(def newtrial (update-ys-pg (update-times (update-G midtrial))))

(def t (keys (:AB (:y newtrial))))
(def y (map first (vals (:AB (:y newtrial)))))
(def omegas (map second (vals (:AB (:y newtrial)))))

(def K (matern-cov-matrix t (:gp-var ABtrials) (:gp-time-scale ABtrials)))
(def obsvars (diagonal-matrix (div 1 omegas)))
(def covmat (add K obsvars (fill (zero-matrix (count t) (count t)) (:gp-mean trial))))
(def meanvector (zero-vector (count t)))
(def loglik1 (d/log-pdf (d/mvnormal meanvector covmat) y))
(def loglik0 (d/log-pdf (d/mvnormal meanvector (add (fill (zero-matrix (count t) (count t)) 1) obsvars)) y))
(def odds (exp (- loglik1 loglik0)))
(/ odds (+ 1 odds ))))))
;(plot/list-plot (:AB (:y newtrial)))
;(switch-probability gp-var gp-time-scale newtrial)



(def my-nested-hashmap {:a "A" :b "B" :c "C" :d "D" :q {:x 1 :y "Y" :z "Z"}})

(defn funzz [{{x :x, y :y} :q, :as state} ] (assoc-in state [:q :x] (inc x) ))
(funzz my-nested-hashmap)





(def reducer (fn [p [m n]] (d/posterior [m] (d/normal n :s2) p)))
(reduce reducer (d/inverse-gamma 3 2) [[5 3] [3 2]])

(def my-maps [{:A 15 :B [14 24] :D {:F 40 :trials ["A" "B"]} } {:A 1 :B [1 2] :D {:F 10 :trials ["C" "D"]} } {:A 2 :B [3 4] :D {:F 20 :trials ["E" "F"]}}])

(defn trial-replace [m]
  (if (map? m)
    (if (contains? m :trials)
      (assoc m :trials (zipmap (map #(keyword (str %)) (range 0 (count (:trials m)))) (:trials m)))
      m)
    m))

(defn -first [x]
  (if (sequential? x) (first x) nil))

(defn furl [x y]
  (if (= (type y) (type (-first x))) (conj x y) (vector x y)))


(defn merge*
  [res lat]
  (if (map? res)
    (merge-with merge* res lat)
    (furl res lat)))

(def foo (apply merge-with merge* (take 3 iterates2)))

(count (first (:G (:0 (:trials (:AB foo))))))

(keyword (str 0))

(type 3.2)
(def iterates2 (apply merge-with merge* (map #(prewalk trial-replace %) (take 900 iterates))))
(:gp-mean (:A iterates2))
(:gp-mean (:A (last iterates)))
(first (iterate inc 1))














(def state (transition initial-state))
(def typ :A)
(def gp-var (:gp-var (:A state)))
(def gp-time-scale (:gp-time-scale (:A state)))
(def single-ys (map :y (get-in state [typ :trials])))
(def dual-ys (map (fn [x] (get-in x [:y typ])) (get-in state [:AB :trials])))
(def y (apply merge (concat single-ys dual-ys)))
(def times (keys y))
(def obs (map first (vals y)))
(def obs-var (map second (vals y)))

(def F (:F (:A (update-single-F times obs obs-var typ state))))
(count (deref (:data (meta F))))
(count y)
(- (d/log-pdf (variance-full-conditional-approx F 0.1) 1.54)
   (d/log-pdf (variance-full-conditional-approx F 0.1) 1.53))
(- (logpdf-prior-gp-approx F 1.54 0.1)
   (logpdf-prior-gp-approx F 1.53 0.1))

(def F (:F (:A new-state)))
(def F (sample-gp 1 gp-time-scale))
(last (map F (range 0 1 0.001)))
(time (logpdf-prior-gp-long F 1 1))
(time (logpdf-prior-gp-test F 1 1))
(time (logpdf-prior-gp-approx F 1 1))
(time (log-likelihood times obs obs-var 1 1))
(keys (:A new-state))
(defn loglikapprox [x] (logpdf-prior-gp-approx F gp-var x))
(defn logliklong [x] (logpdf-prior-gp-long F gp-var x))
(defn loglikmar [x] (log-likelihood times obs obs-var gp-var x))
(time (normalize-log (map loglikapprox [0.1 0.2 0.3 0.4 0.5])))

(defn normalize-log [coll]
  (let [a (reduce max coll)
        shifted (map - coll (repeat a))
        expshifted (map exp shifted)
        Z (reduce + expshifted)
        ]
    (map #(/ % Z) expshifted)))

(def dist (variance-full-conditional-approx F gp-time-scale))
(ic/view (scatter-plot (range 0.1 2 0.1) (normalize-log (map (d/log-pdf dist) (range  0.1 2 0.1) ))))
(ic/view (scatter-plot (range 0.1 2 0.1) (normalize-log (map #(logpdf-prior-gp-approx F % gp-time-scale) (range  0.1 2 0.1) ))))
(ic/view (scatter-plot (range 0.1 2 0.1) (normalize-log (map #(logpdf-prior-gp-long F % gp-time-scale) (range  0.1 2 0.1) ))))
(ic/view (scatter-plot (range 0.01 0.3 0.01) (normalize-log (map loglikapprox (range  0.01 0.3 0.01) ))))
(ic/view (scatter-plot (range 0.01 0.3 0.01) (normalize-log (map logliklong (range  0.01 0.3 0.01) ))))
(ic/view (scatter-plot (range 0.01 0.3 0.01) (normalize-log (map loglikmar (range  0.01 0.3 0.01) ))))
(loglik 1.5)
(loglik 0.4)
(- (logpdf-mvnormal F AG Q) (logpdf-mvnormal F AG Q) )
(logpdf-normal (mget F 2) (mget AG 2) (mget Q 2 2))
(logpdf-mvnormal F AG Q)
(logpdf-normal (mget F 2) (mget AG 2) (mget Q 2 2))
