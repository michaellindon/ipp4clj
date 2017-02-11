(ns ipp4clj.multiplexing
  (:require [clojure.core.matrix :refer :all]
            [incanter.core :as ic]
            [incanter.stats :refer :all]
            [incanter.charts :as ip]
            [clojure.data.avl :as avl]
            [ssm4clj.core :refer :all]
            [ssm4clj.misc :refer :all]
            [ipp4clj.misc :refer :all]))

(defn zero-fn [x] 0.0)
(defn group-bool [pred coll]
  (let [bool-map (group-by pred coll)]
    [(get bool-map true) (get bool-map false)]))

(defn generate-single-trials [num-trials gp-mean gp-var gp-time-scale max-intensity]
  (let [F (sample-gp {} gp-var gp-time-scale)
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean ) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))
        obs-times (sample-ppp intensity max-intensity 0 1)
        aug-times (sample-ppp thin-intensity max-intensity 0 1)
        obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
        aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
        obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
        aug-map (into (avl/sorted-map) (map vector aug-times aug-y))
        ]
    {:F F :gp-mean gp-mean :gp-var gp-var :gp-time-scale gp-time-scale :max-intensity max-intensity
     :trials (take num-trials (repeatedly (fn [] {:obs-times obs-times
                                                  :aug-times aug-times
                                                  :y (merge obs-map aug-map)
                                                  :start-t 0
                                                  :end-t 1})))}))

(defn generate-dual-trials [num-trials gp-var gp-time-scale switching-prob A B]
  (let [{F-A :F gp-mean-A :gp-mean max-intensity :max-intensity} A
        {F-B :F gp-mean-B :gp-mean max-intensity :max-intensity} B
        f-A (comp first F-A)
        m-A (comp (partial + gp-mean-A) f-A)
        intensity-A (comp (partial * max-intensity) probit m-A)
        f-B (comp first F-B)
        m-B (comp (partial + gp-mean-B ) f-B)
        intensity-B (comp (partial * max-intensity) probit m-B)]
    {:gp-var gp-var :gp-time-scale gp-time-scale
     :trials (take
              num-trials
              (repeatedly
               (fn [] (let [G (sample-gp {} gp-var gp-time-scale)
                            switching (sample-binomial 1 :size 1 :prob switching-prob)
                            g (if (= switching 1) (comp first G) zero-fn)
                            gp-mean-AB (sample-normal 1)
                            alpha (comp probit (partial + gp-mean-AB) g)
                            intensity (fn [t] (+ (* (alpha t) (intensity-A t))
                                                 (* (- 1 (alpha t)) (intensity-B t))))
                            obs-times (sample-ppp intensity max-intensity 0 1)
                            intensity-Aar (fn [x] (* (- 1 (alpha x)) (intensity-A x)))
                            Aar (sample-ppp intensity-Aar max-intensity 0 1)
                            intensity-Arr (fn [x] (- max-intensity (intensity-A x)))
                            Arr (sample-ppp intensity-Arr max-intensity 0 1)
                            intensity-Bar (fn [x] (* (alpha x) (intensity-B x)))
                            Bar (sample-ppp intensity-Bar max-intensity 0 1)
                            intensity-Brr (fn [x] (- max-intensity (intensity-B x )))
                            Brr (sample-ppp intensity-Brr max-intensity 0 1)
                            ass-prob (fn [x] (/ (* (- 1 (alpha x)) (intensity-B x))
                                                (+ (* (alpha x) (intensity-A x))
                                                   (* (- 1 (alpha x)) (intensity-B x)))))
                            ass-pred (fn [t] (= 1 (sample-binomial 1 :size 1 :prob (ass-prob t))))
                            [Baa Aaa] (group-bool ass-pred obs-times)
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
                         :switching switching
                         :y {:A y-A :B y-B :AB y-AB}
                         :G G
                         :start-t 0
                         :end-t 1}))))}))

(defn plot-single [x]
  (let [{F :F
         gp-mean :gp-mean
         max-intensity :max-intensity
         trials :trials} x
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean ) f)]
    (do (ic/view (ip/function-plot f 0 1 :title "Gaussian Process f"))
        (ic/view (ip/histogram (apply concat (map :obs-times trials)) :title "PSTH" :nbins 20))
        (ic/view (ip/function-plot intensity 0 1 :title "Intensity")))))

(def Atrials (generate-single-trials 15 0 2 0.2 100))
(def Btrials (generate-single-trials 15 0 2 0.2 100))
(def ABtrials (generate-dual-trials 15 2 0.2 0.5 Atrials Btrials))
(plot-single Btrials)

(def initial-state {:A Atrials :B Btrials :AB ABtrials})

(defn ulogpdf-gamma [a b x]
 (if (neg? x)
     Double/NEGATIVE_INFINITY
     (- (* (dec a) (log x)) (* b x))))
(defn logprior-gp-time-scale [x] (ulogpdf-gamma 2 1 x))
(defn logprior-gp-var [x] (ulogpdf-gamma 2 1 x))

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
    (assoc trial :y {:A y-A :B y-B :AB y-AB})))

(defn update-dual-times [intensity-A intensity-B max-intensity trial]
  (let [{start-t :start-t
         end-t :end-t
         G :G
         gp-mean :gp-mean
         obs-times :obs-times} trial
        alpha (comp probit (partial + gp-mean) first G)
        intensity-Aar (fn [x] (* (- 1 (alpha x)) (intensity-A x)))
        Aar (sample-ppp intensity-Aar max-intensity start-t end-t)
        intensity-Arr (fn [x] (- max-intensity (intensity-A x)))
        Arr (sample-ppp intensity-Arr max-intensity start-t end-t)
        intensity-Bar (fn [x] (* (alpha x) (intensity-B x)))
        Bar (sample-ppp intensity-Bar max-intensity start-t end-t)
        intensity-Brr (fn [x] (- max-intensity (intensity-B x )))
        Brr (sample-ppp intensity-Brr max-intensity start-t end-t)
        ass-prob (fn [x] (/ (* (- 1 (alpha x)) (intensity-B x))
                            (+ (* (alpha x) (intensity-A x))
                               (* (- 1 (alpha x)) (intensity-B x)))))
        ass-pred (fn [t] (= 1 (sample-binomial 1 :size 1 :prob (ass-prob t))))
        [Baa Aaa] (group-bool ass-pred obs-times)]
    (assoc trial :Aar Aar :Arr Arr :Bar Bar :Brr Brr :Baa Baa :Aaa Aaa)))

(defn update-dual-trials [state]
  (let [{F-A :F gp-mean-A :gp-mean max-intensity :max-intensity} (:A state)
        {F-B :F gp-mean-B :gp-mean max-intensity :max-intensity} (:B state)
        f-A (comp first F-A)
        m-A (comp (partial + gp-mean-A ) f-A)
        intensity-A (comp (partial * max-intensity) probit m-A)
        f-B (comp first F-B)
        m-B (comp (partial + gp-mean-B ) f-B)
        intensity-B (comp (partial * max-intensity) probit m-B)
        {trials :trials} (:AB state)
        update-times (partial update-dual-times intensity-A intensity-B max-intensity)
        update-ys (partial update-dual-ys m-A m-B)
        update-trial (comp update-ys update-times)
        ]
    (assoc-in state [:AB :trials] (map update-trial trials))))

(defn dual-log-likelihood [gp-time-scale gp-var trial]
  (let [y (get-in trial [:y :AB])
        times (keys y)
        obs (vals y)]
    (log-likelihood times obs 1 gp-var gp-time-scale)))

(defn update-dual-gp-time-scale [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)
        conditional (fn [x]
                      (+ (logprior-gp-time-scale x)
                         (reduce + 0
                                 (map (partial dual-log-likelihood x gp-var) trials))))]
    (assoc-in state [:AB :gp-time-scale] (sample-slice conditional 1 gp-time-scale))))

(defn update-dual-gp-var [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)
        conditional (fn [x]
                      (+ (logprior-gp-var x)
                         (reduce + 0
                                 (map (partial dual-log-likelihood gp-time-scale x) trials))))]
    (assoc-in state [:AB :gp-var] (sample-slice conditional 1 gp-var))))


(defn dual-gp-mean [gp-var gp-time-scale trial]
  (let [{y :AB} (:y trial)
        times (keys y)
        obs (vals y)
        mus2 (mean-conditional times obs 1 gp-var gp-time-scale)
        new-mean (sample-normal 1 :mean (first mus2) :sd (sqrt (second mus2)))]
    (assoc trial :gp-mean new-mean)))

(defn update-dual-gp-mean [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)]
    (assoc-in state [:AB :trials] (map (partial dual-gp-mean gp-var gp-time-scale) trials))))

(defn dual-G [gp-var gp-time-scale trial]
  (let [{y :AB} (:y trial)
        {gp-mean :gp-mean} trial
        times (keys y)
        obs (vals y)
        centered-obs (map (fn [x] (- x gp-mean)) obs)
        pivots (FFBS times centered-obs 1 gp-var gp-time-scale)]
    (assoc trial :G (sample-gp pivots gp-var gp-time-scale))))

(defn update-dual-G [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)
        foo (map (partial dual-G gp-var gp-time-scale) (:trials (:AB state)))
        ]
    (assoc-in state [:AB :trials] foo)))

(defn update-single-ys [gp-mean f trial]
  (let [{obs-times :obs-times
         aug-times :aug-times} trial
        obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
        aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
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
        intensity (comp (partial * max-intensity) probit (partial + gp-mean) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))
        update-ys (partial update-single-ys gp-mean f)
        update-single-times (partial update-single-times thin-intensity max-intensity)
        update-trial (comp update-ys update-single-times)
        ]
    (assoc-in state [typ :trials] (map update-trial trials))))

(defn update-single-gp-mean [times obs typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale} (typ state)
        mus2 (mean-conditional times obs 1 gp-var gp-time-scale)
        new-mean (sample-normal 1 :mean (first mus2) :sd (sqrt (second mus2)))]
    (assoc-in state [typ :gp-mean] new-mean)))

(defn update-single-gp-time-scale [times obs typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale} (typ state)
        conditional (fn [x] (+ (logprior-gp-time-scale x)
                               (log-likelihood times obs 1 gp-var x)))]
    (assoc-in state [typ :gp-time-scale] (sample-slice conditional 1 gp-time-scale))))

(defn update-single-gp-var [times obs typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale} (typ state)
        conditional (fn [x] (+ (logprior-gp-var x)
                               (log-likelihood times obs 1 x gp-time-scale)))]
    (assoc-in state [typ :gp-var] (sample-slice conditional 1 gp-var))))

(defn update-single-F [times obs typ state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         gp-mean :gp-mean} (typ state)
        centered-obs (map (fn [x] (- x gp-mean)) obs)
        pivots (FFBS times centered-obs 1 gp-var gp-time-scale)]
    (assoc-in state [typ :F] (sample-gp pivots gp-var gp-time-scale))))

(defn update-single-intensity [typ state]
  (let [single-ys (map :y (get-in state [typ :trials]))
        dual-ys (map (fn [x] (get-in x [:y typ])) (get-in state [:AB :trials]))
        y (apply merge (concat single-ys dual-ys))
        times (keys y)
        obs (vals y)
        update-gp-mean (partial update-single-gp-mean times obs typ)
        update-gp-time-scale (partial update-single-gp-time-scale times obs typ)
        update-gp-var (partial update-single-gp-var times obs typ)
        update-F (partial update-single-F times obs typ)
        update (comp update-F
                     update-gp-var
                     update-gp-time-scale
                     update-gp-mean)]
    (update state)))


(defn test-update-dual-trials [state]
  (let [{F-A :F gp-mean-A :gp-mean max-intensity :max-intensity} (:A state)
        {F-B :F gp-mean-B :gp-mean max-intensity :max-intensity} (:B state)
        f-A (comp first F-A)
        m-A (comp (partial + gp-mean-A ) f-A)
        intensity-A (comp (partial * max-intensity) probit m-A)
        f-B (comp first F-B)
        m-B (comp (partial + gp-mean-B ) f-B)
        intensity-B (comp (partial * max-intensity) probit m-B)
        {trials :trials} (:AB state)
        update-ys (partial update-dual-ys m-A m-B)
        update-trial update-ys
        ]
    (assoc-in state [:AB :trials] (map update-trial trials))))

(def transition (comp
                 update-dual-trials
                 (partial update-single-intensity :A)
                 (partial update-single-trials :A)
                 (partial update-single-intensity :B)
                 (partial update-single-trials :B)))


(def mcmc (iterate transition initial-state))
(def iterates (doall (take 200 mcmc)))

(defn plot-intensity [typ iterates truth]
  (let[
       mean-function (fn [x] (mean (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates))))
       lower-function (fn [x] (quantile (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates)) :probs 0.025))
       upper-function (fn [x] (quantile (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates)) :probs 0.975))
       ribbon (ip/function-plot mean-function 0 1)
       ]
    (do
      (.setSeriesPaint (.getRenderer (.getPlot ribbon) 0) 0 java.awt.Color/red)
      (ip/add-function ribbon upper-function 0 1 :series-label 1)
      (.setSeriesPaint (.getRenderer (.getPlot ribbon) 1) 0 java.awt.Color/pink)
      (ip/add-function ribbon lower-function 0 1 :series-label 2)
      (.setSeriesPaint (.getRenderer (.getPlot ribbon) 2) 0 java.awt.Color/pink)
      (ip/add-function ribbon (fn [x] (* (:max-intensity truth) (probit (+ (:gp-mean truth) ((comp first (:F truth)) x))))) 0 1 :series-label 3)
      (.setSeriesPaint (.getRenderer (.getPlot ribbon) 3) 0 java.awt.Color/black)
      (ic/view ribbon)
      )))

(plot-intensity :A iterates Atrials)

(defn mean-function [x] (mean (map (fn [y] (+ (:gp-mean y) (first ((:F y) x)))) (map :B iterates))))
(defn lower-function [x] (quantile (map (fn [y] (+ (:gp-mean y) (first ((:F y) x)))) (map :B iterates)) :probs 0.025))
(defn upper-function [x] (quantile (map (fn [y] (+ (:gp-mean y) (first ((:F y) x)))) (map :B iterates)) :probs 0.975))
(def ribbon (ip/function-plot mean-function 0 1))
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 0) 0 java.awt.Color/red)
(ip/add-function ribbon upper-function 0 1 :series-label 1)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 1) 0 java.awt.Color/pink)
(ip/add-function ribbon lower-function 0 1 :series-label 2)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 2) 0 java.awt.Color/pink)
(ip/add-function ribbon (fn [x] (+ (:gp-mean Btrials) ((comp first (:F Btrials)) x))) 0 1 :series-label 3)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 3) 0 java.awt.Color/black)
(ic/view ribbon)

(ic/view (ip/time-series-plot (range 0 (ic/length iterates)) (map :gp-mean iterates)))
(ic/view (ip/histogram (map :gp-mean iterates) :nbins 100))
(ic/view (ip/histogram (map :gp-var iterates) :nbins 100))
(ic/view (ip/histogram (map :gp-time-scale iterates) :nbins 100))

(map (partial dual-G 1 1) (:trials (:AB initial-state)))
((partial dual-G 1 1) (first (:AB initial-state)))
(map (partial dual-G 1 1) (:AB initial-state))


(keys (:AB initial-state))
(update-single-intensity :A (update-single-trials :A initial-state))
partition-by
filter
(group-by #(> 4 %) (range 1 10))
(get  (group-by #(> 4 %) (range 1 10)) true)
(group-bool #(> 4 %) (range 1 10))

(assoc {:foo 11} :a 3 :b 5 :d 3)

(concat [[99 100] [100 111]] [[1 1] [2 3] [4 5]] [[14 5] [5 6] [7 8]])

(reduce (fn [m x] (assoc m x (inc x))) (avl/sorted-map) [1 2 3])

(reduce + 0 [])
(mean-function 0.3)

(map (partial dual-G 2 0.2) (:trials (:AB initial-state)))

(dual-G 2 0.2 (first (:trials (:AB initial-state))))

(map (partial dual-G 1 1) (:trials (:AB initial-state)))

initial-state
(update-dual-G initial-state)
(def updated-state (update-dual-G initial-state))
((:G  (dual-G 1 1 (first (:trials (:AB initial-state))))) 0)
(:trials (:AB updated-state))

(:gp-mean (first (:trials (:AB initial-state))))
(:gp-mean (first (:trials (:AB updated-state))))

( (:G (first (:trials (:AB initial-state)))) 0)
( (:G (first (:trials (:AB updated-state)))) 0)
