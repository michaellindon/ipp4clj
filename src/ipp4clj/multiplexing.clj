(ns ipp4clj.multiplexing
  (:require [clojure.core.matrix :refer :all]
            [incanter.stats :refer :all]
            [clojure.data.avl :as avl]
            [ssm4clj.core :refer :all]
            [ssm4clj.misc :refer :all]
            [ipp4clj.misc :refer :all]
            [gorilla-plot.core :as plot]))

(defn zero-fn [x] 0.0)
(defn zero-3fn [x] (matrix [0.0 0.0 0.0]))
(defn group-bool [pred coll]
  (let [bool-map (group-by pred coll)]
    [(get bool-map true) (get bool-map false)]))

(defn generate-single-trials [num-trials gp-mean gp-var gp-time-scale max-intensity]
  (let [F (sample-gp {} gp-var gp-time-scale)
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean ) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))]
    {:F F :gp-mean gp-mean :gp-var gp-var :gp-time-scale gp-time-scale :max-intensity max-intensity
     :trials
     (take
      num-trials
      (repeatedly (fn []
                    (let [obs-times (sample-ppp intensity max-intensity 0 1)
                          aug-times (sample-ppp thin-intensity max-intensity 0 1)
                          obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
                          aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
                          obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
                          aug-map (into (avl/sorted-map) (map vector aug-times aug-y))]
                        {:obs-times obs-times
                         :aug-times aug-times
                         :y (merge obs-map aug-map)
                         :start-t 0
                         :end-t 1}))))}))

(defn intensity-functions
  "Provides all relevant intensity functions for AB trials.
  A modular function useful for unit testing"
  [alpha intensity-A intensity-B max-intensity]
  (let [i-Aar (fn [x] (* (- 1 (alpha x)) (intensity-A x)))
        i-Arr (fn [x] (- max-intensity (intensity-A x)))
        i-Bar (fn [x] (* (alpha x) (intensity-B x)))
        i-Brr (fn [x] (- max-intensity (intensity-B x )))
        i-Aaa (fn [x] (* (alpha x) (intensity-A x)))
        i-Baa (fn [x] (* (- 1 (alpha x)) (intensity-B x)))]
    [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr]))

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
               (fn [] (let [
                            switching (sample-binomial 1 :size 1 :prob switching-prob)
                            G (if (= switching 1) (sample-gp {} gp-var gp-time-scale) zero-3fn)
                            g (comp first G)
                            gp-mean-AB (sample-normal 1)
                            alpha (comp probit (partial + gp-mean-AB) g)
                            [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr] (intensity-functions
                                                                   alpha
                                                                   intensity-A
                                                                   intensity-B
                                                                   max-intensity)
                            start-t 0
                            end-t 1
                            Aaa (sample-ppp i-Aaa max-intensity start-t end-t)
                            Aar (sample-ppp i-Aar max-intensity start-t end-t)
                            Arr (sample-ppp i-Arr max-intensity start-t end-t)
                            Baa (sample-ppp i-Baa max-intensity start-t end-t)
                            Bar (sample-ppp i-Bar max-intensity start-t end-t)
                            Brr (sample-ppp i-Brr max-intensity start-t end-t)
                            obs-times (concat Aaa Baa)
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
                         :start-t start-t
                         :end-t end-t}))))}))

(defn ulogpdf-gamma [a b x]
 (if (neg? x)
     Double/NEGATIVE_INFINITY
     (- (* (dec a) (log x)) (* b x))))
(defn logprior-gp-time-scale [x] (ulogpdf-gamma 2 5 x))
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
        [i-Aaa i-Aar i-Arr i-Baa i-Bar i-Brr] (intensity-functions
                                               alpha
                                               intensity-A
                                               intensity-B
                                               max-intensity)
        Aar (sample-ppp i-Aar max-intensity start-t end-t)
        Arr (sample-ppp i-Arr max-intensity start-t end-t)
        Bar (sample-ppp i-Bar max-intensity start-t end-t)
        Brr (sample-ppp i-Brr max-intensity start-t end-t)
        AB-prob (fn [x] (/ (i-Aaa x) (+ (i-Aaa x) (i-Baa x))))
        AB-pred (fn [t] (= 1 (sample-binomial 1 :size 1 :prob (AB-prob t))))
        [Aaa Baa] (group-bool AB-pred obs-times)]
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

(defn dual-log-likelihood
  [gp-time-scale gp-var trial]
  (let [y (get-in trial [:y :AB])
        times (keys y)
        obs (vals y)]
    (log-likelihood times obs 1 gp-var gp-time-scale)))

(defn switch-log-likelihood
  [gp-time-scale gp-var trial]
 (let [y (get-in trial [:y :AB])
       times (keys y)
       obs (vals y)
       {switching :switching
        gp-mean :gp-mean} trial]
   (if (= switching 1)
     (log-likelihood times obs 1 gp-var gp-time-scale)
     (reduce + 0 (map logpdf-normal obs (repeat gp-mean) (repeat 1.0))))))

(dual-log-likelihood 1.0 1.0 (first (:trials ABtrials)))
(dual-log-likelihood 0 0.2 1.0 1.0 (first (:trials ABtrials)))

(defn dual-switching )


(defn update-dual-gp-time-scale [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)
        conditional (fn [x]
                      (+ (logprior-gp-time-scale x)
                         (reduce + 0
                                 (map (partial dual-log-likelihood x gp-var) trials))))]
    (assoc-in state [:AB :gp-time-scale] (sample-slice conditional 1 gp-time-scale))))

;;TODO SHOULD THIS BE MEAN CENTERED
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
  (if (= (:switching trial)) 1
      (let [{y :AB} (:y trial)
            times (keys y)
            obs (vals y)
            mus2 (mean-conditional times obs 1 gp-var gp-time-scale)
            new-mean (sample-normal 1 :mean (first mus2) :sd (sqrt (second mus2)))]
        (assoc trial :gp-mean new-mean))
      ;;TODO IMPLEMENT DRAW FROM UNIVARIATE NORMAL
      ))

(defn update-dual-gp-mean [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)]
    (assoc-in state [:AB :trials] (map (partial dual-gp-mean gp-var gp-time-scale) trials))))

(defn dual-G [gp-var gp-time-scale trial]
  (let [{switching :switching} trial]
    (if (= switching 0)
      (assoc trial :G zero-3fn)
      (let [{y :AB} (:y trial)
            {gp-mean :gp-mean} trial
            times (keys y)
            obs (vals y)
            centered-obs (map (fn [x] (- x gp-mean)) obs)
            pivots (FFBS times centered-obs 1 gp-var gp-time-scale)]
        (assoc trial :G (sample-gp pivots gp-var gp-time-scale))))))

(defn update-dual-G [state]
  (let [{gp-var :gp-var
         gp-time-scale :gp-time-scale
         trials :trials} (:AB state)]
    (assoc-in state [:AB :trials] (map (partial dual-G gp-var gp-time-scale) trials))))

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

(def transition (comp
                 update-dual-G
                 update-dual-gp-time-scale
                 update-dual-gp-var
                 update-dual-trials
                 (partial update-single-intensity :B)
                 (partial update-single-trials :B)
                 (partial update-single-intensity :A)
                 (partial update-single-trials :A)))

(def Atrials (generate-single-trials 6 0 2 0.2 100))
(def Btrials (generate-single-trials 6 0 2 0.2 100))
(def ABtrials (generate-dual-trials 6 2 0.2 0.5 Atrials Btrials))

(def initial-state {:A Atrials :B Btrials :AB ABtrials})
(def mcmc (iterate transition initial-state))
(def iterates (doall (take 100 (drop 10 mcmc))))

(defn plot-intensity
  ([typ iterates]
   (let[
        ;mean-function (fn [x] (mean (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates))))
        lower-function (fn [x] (quantile (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates)) :probs 0.025))
        upper-function (fn [x] (quantile (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates)) :probs 0.975))
        ]
     (plot/compose (plot/plot upper-function [0 1] :plot-range [[0 1] [0 100]] :colour "skyblue")
                   (plot/plot lower-function [0 1] :plot-range [[0 1] [0 100]] :colour "skyblue"))
     ))
  ([typ iterates truth]
   (let [true-function (fn [x] (* (:max-intensity truth) (probit (+ (:gp-mean truth) ((comp first (:F truth)) x)))))]
     (plot/compose
      (plot-intensity typ iterates)
      (plot/plot true-function [0 1] :plot-range [[0 1] [0 100]] :colour "#FA8072"))
     )
   ))

(defn plot-gp
  ([typ iterates]
   (let[
        ;mean-function (fn [x] (mean (map (fn [y] (* (:max-intensity y) (probit (+ (:gp-mean y) (first ((:F y) x)))))) (map typ iterates))))
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
  (let [max-intensity (map #(get-in % [:A :max-intensity]) iterates)
        f-A (map #(comp first (get-in % [:A :F])) iterates)
        gp-mean-A (map #(get-in % [:A :gp-mean]) iterates)
        f-B (map #(comp first (get-in % [:B :F])) iterates)
        gp-mean-B (map #(get-in % [:B :gp-mean]) iterates)
        intensity-A (map (fn [Lambda gp-mean f] (fn [t] (* Lambda (probit (+ gp-mean (f t)))))) max-intensity gp-mean-A f-A)
        intensity-B (map (fn [Lambda gp-mean f] (fn [t] (* Lambda (probit (+ gp-mean (f t)))))) max-intensity gp-mean-B f-B)]
    {:intensity-A intensity-A :intensity-B intensity-B}
    ))

(defn plot-dual [i iterates truth]
  (let [g-AB (map #(comp first (:G (nth (get-in % [:AB :trials]) i))) iterates)
        gp-mean-AB (map #(:gp-mean (nth (get-in % [:AB :trials]) i)) iterates)
        mean-histogram (plot/histogram gp-mean-AB :bins 30)
        mean-traceplot (plot/list-plot gp-mean-AB :joined true)
        g-plot (plot-gp-dual g-AB i truth)
        intensity-plot (plot-intensity-dual i iterates truth gp-mean-AB g-AB)
        ]
    [mean-histogram mean-traceplot g-plot intensity-plot]))

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

(defn plot-intensity-dual
  ([i iterates gp-mean-AB g-AB]
   (let [{intensity-A :intensity-A
          intensity-B :intensity-B} (extract-params iterates)
         alpha (map (fn [gp-mean g] (fn [t] (probit (+ gp-mean (g t))))) gp-mean-AB g-AB)
         intensity-AB (map (fn [in-A in-B a]
                             (fn [t] (+ (* (a t) (in-A t))
                                        (* (- 1 (a t)) (in-B t))))) intensity-A intensity-B alpha)
         [lower-intensity upper-intensity] (quartile-functions intensity-AB)]
     (plot/compose (plot/plot upper-intensity [0 1] :plot-range [[0 1] [0 100]] :colour "skyblue")
                   (plot/plot lower-intensity [0 1] :colour "skyblue"))
     ))
  ([i iterates truth gp-mean-AB g-AB]
   (let [{intensity-AB :intensity-AB} (get-truth i truth)]
     (plot/compose
      (plot/plot intensity-AB [0 1] :plot-range [[0 1] [0 100]] :colour "#FA8072")
      (plot-intensity-dual i iterates gp-mean-AB g-AB)))
   ))


(defn get-truth [i truth]
  (let [max-intensity (get-in truth [:A :max-intensity])
        f-A (comp first (get-in truth [:A :F]))
        gp-mean-A (get-in truth [:A :gp-mean])
        f-B (comp first (get-in truth [:B :F]))
        gp-mean-B (get-in truth [:B :gp-mean])
        g-AB (comp first (:G (nth (get-in truth [:AB :trials]) i)))
        gp-mean-AB (:gp-mean (nth (get-in truth [:AB :trials]) i))
        intensity-A (fn [t] (* max-intensity (probit (+ gp-mean-A (f-A t)))))
        intensity-B (fn [t] (* max-intensity (probit (+ gp-mean-B (f-B t)))))
        alpha (fn [t] (probit (+ gp-mean-AB (g-AB t))))
        intensity-AB (fn [t] (+ (* (alpha t) (intensity-A t))
                                (* (- 1 (alpha t)) (intensity-B t))))]
    {:intensity-AB intensity-AB
     :g-AB g-AB}
    ))


(defn plot-gp-duals [iterates true-Gs]
  (let [Gstreams (for [i (range 0 (count (:trials ABtrials)))]
                   (map #(:G (nth (get-in % [:AB :trials]) i)) iterates))]
    (map plot-gp-dual Gstreams true-Gs)))

(plot-gp :A iterates Atrials)

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
   [(latex-view "\\tau_{AB}") 
    (str "Characteristic time-scale of all Gaussian Processes g_AB. Truth " (:gp-time-scale ABtrials)) 
    (plot-param :gp-time-scale :AB iterates)
    (plot-trace :gp-time-scale :AB iterates)]
   [(latex-view "\\sigma^2_{AB}") 
    (str "Variance of all Gaussian Processes g_AB. Truth " (:gp-var ABtrials)) 
    (plot-param :gp-var :AB iterates)
    (plot-trace :gp-var :AB iterates)] ]
    :columns 
  ["Paramter" "Description" "Posterior Summary" "TracePlot"]) )
