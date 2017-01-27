(use '(incanter stats))
(use 'clojure.core.matrix)
(require '[incanter.core :as ic])
(require '[incanter.charts :as ip])
(require '[clojure.data.avl :as avl])
(require '[clojure.core.matrix.linear :as la])
(require '[clojure.math.combinatorics :as combo])
(require '[ipp4clj.statespace :refer :all])
(require '[ipp4clj.misc :refer :all])
(require '[ipp4clj.gaussianprocess :refer :all])

(defn generate-single-trials [num-trials gp-mean gp-var gp-time-scale max-intensity]
  (let [F (sample-gp {} gp-var gp-time-scale)
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean ) f)]
    {:F F :gp-mean gp-mean :gp-var gp-var :gp-time-scale gp-time-scale :max-intensity max-intensity
     :trials (take num-trials (repeatedly (fn [] {:obs-times (sample-ppp intensity max-intensity 0 1)
                                                  :aug-times []
                                                  :y (avl/sorted-map)
                                                  :start-t 0
                                                  :end-t 1})))}))

(defn zero-fn [x] 0.0)

(defn generate-dual-trials [num-trials gp-var gp-time-scale switching-prob A B]
  (let [{F-A :F gp-mean-A :gp-mean max-intensity :max-intensity} A
        {F-B :F gp-mean-B :gp-mean max-intensity :max-intensity} B
        f-A (comp first F-A)
        intensity-A (comp (partial * max-intensity) probit (partial + gp-mean-A ) f-A)
        f-B (comp first F-B)
        intensity-B (comp (partial * max-intensity) probit (partial + gp-mean-B ) f-B)]
    {:gp-var gp-var :gp-time-scale gp-time-scale
     :trials (take
              num-trials
              (repeatedly (fn [] (let [G (sample-gp {} gp-var gp-time-scale)
                                       switching (sample-binomial 1 :size 1 :prob switching-prob)
                                       g (if (= switching 1) (comp first G) zero-fn)
                                       gp-mean (sample-normal 1)
                                       alpha (comp probit (partial + gp-mean) g)
                                       intensity (fn [t] (+ (* (alpha t) (intensity-A t))
                                                            (* (- 1 (alpha t)) (intensity-B t))))]
                                   {:obs-times (sample-ppp intensity max-intensity 0 1)
                                    :aug-times []
                                    :switching switching
                                    :y (avl/sorted-map)
                                    :G G
                                    :start-t 0
                                    :end-t 1}))))}))

(defn plot-single [x]
  (let [{F :F
         gp-mean :gp-mean
         max-intensity :max-intensity
         trials :trials} x
        {obs-times :obs-times} trials
        spike-train-plot (ip/scatter-plot (first obs-times) (replicate (count (first obs-times)) 0))
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean ) f)]
    (do (ic/view (ip/function-plot f 0 1 :title "Gaussian Process f"))
        (ic/view spike-train-plot)
        (ic/view (ip/function-plot intensity 0 1 :title "Intensity")))))

(def spike-train-plot (ip/scatter-plot (sample-normal 10) (replicate 10 1)))
(ip/add-points spike-train-plot (sample-normal 10) (replicate 10 2))
(ic/view spike-train-plot)

(def Atrials (generate-single-trials 3 0 2 0.2 1000))
(def Btrials (generate-single-trials 3 0 2 0.2 1000))
(def ABtrials (generate-dual-trials 3 2 0.2 0.5 Atrials Btrials))
(plot-single Atrials)

(def initial-state {:A {:F F
                        :gp-var gp-var
                        :gp-time-scale gp-time-scale
                        :gp-mean gp-mean
                        :max-intensity max-intensity
                        :obs-times obs
                        :aug-times aug
                        :y (merge obs-map aug-map)
                        :start-t 0
                        :end-t 1}})

(defn ulogpdf-gamma [a b x]
 (if (neg? x)
     Double/NEGATIVE_INFINITY
     (- (* (dec a) (log x)) (* b x))))
(defn logprior-gp-time-scale [x] (ulogpdf-gamma 2 1 x))
(defn logprior-gp-var [x] (ulogpdf-gamma 2 1 x))

(first (:obs-times (first Atrials)))
(first (:obs-times (second Atrials)))
(type Atrials)

(defn update-single-trials [typ state]
  (let [{F :F
         max-intensity :max-intensity
         gp-mean :gp-mean
         trials :trials} (typ state)
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))
        update-ys (fn [trial]
                   (let [{obs-times :obs-times
                          aug-times :aug-times} trial
                         obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
                         aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)]
                     (assoc trial :y (merge obs-map aug-map))))
        update-aug-times (fn [trial]
                          (let [{start-t :start-t
                                 end-t :end-t} trial]
                            (assoc trial :aug-times (sample-ppp thin-intensity max-intensity start-t end-t))))
        update-trial (comp update-ys update-aug-times)
        ]
    (assoc-in state [typ :trials] (map update-trial trials))))

(defn update-single-intensity [typ state]
  (let [single-ys (map :y (get-in state [typ :trials]))
        dual-ys (map :y (get-in state [:AB :trials]))
        y (apply merge (concat single-ys dual-ys))
        times (keys y)
        obs (vals y)
        update-gp-mean (fn [typ state]
                         (let [{gp-var :gp-var
                                gp-time-scale :gp-time-scale} (typ state)
                               mus2 (mean-conditional times obs 1 gp-var gp-time-scale)
                               new-mean (sample-normal 1 :mean (first mus2) :sd (sqrt (second mus2)))]
                           (assoc-in state [typ :gp-mean] new-mean)))
        update-gp-time-scale (fn [typ state]
                               (let [{gp-var :gp-var
                                      gp-time-scale :gp-time-scale} (typ state)
                                     conditional (fn [x] (+ (logprior-gp-time-scale x)
                                                            (log-likelihood times obs 1 gp-var x)))]
                                 (assoc-in state [typ :gp-time-scale] (sample-slice conditional 1 gp-time-scale))))
        update-gp-var (fn [typ state]
                        (let [{gp-var :gp-var
                               gp-time-scale :gp-time-scale} (typ state)
                              conditional (fn [x] (+ (logprior-gp-var x)
                                                     (log-likelihood times obs 1 x gp-time-scale)))]
                          (assoc-in state [typ :gp-var] (sample-slice conditional 1 gp-var))))
        update-F (fn [typ state]
                   (let [{gp-var :gp-var
                          gp-time-scale :gp-time-scale
                          gp-mean :gp-mean} (typ state)
                         centered-obs (map (fn [x] (- x gp-mean)) obs)
                         pivots (FFBS times centered-obs 1 gp-var gp-time-scale)]
                     (assoc-in state [typ :F] (sample-gp pivots gp-var gp-time-scale))))
        update (comp (partial update-F typ)
                     (partial update-gp-var typ)
                     (partial update-gp-time-scale typ)
                     (partial update-gp-mean typ))]
    (update state)))

(ic/view (ip/plot (sample-normal 10) ones ))


(def mcmc (iterate transition initial-state))
(def iterates (take 200 mcmc))
(defn mean-function [x] (mean (map (fn [y] (+ (:gp-mean y) (first ((:F y) x)))) iterates)))
(defn lower-function [x] (quantile (map (fn [y] (+ (:gp-mean y) (first ((:F y) x)))) iterates) :probs 0.025))
(defn upper-function [x] (quantile (map (fn [y] (+ (:gp-mean y) (first ((:F y) x)))) iterates) :probs 0.975))
(def ribbon (ip/function-plot mean-function 0 1))
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 0) 0 java.awt.Color/red)
(ip/add-function ribbon upper-function 0 1 :series-label 1)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 1) 0 java.awt.Color/pink)
(ip/add-function ribbon lower-function 0 1 :series-label 2)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 2) 0 java.awt.Color/pink)
(ip/add-function ribbon (fn [x] (+ gp-mean (f x))) 0 1 :series-label 3)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 3) 0 java.awt.Color/black)
(ic/view ribbon)

(ic/view (ip/time-series-plot (range 0 (ic/length iterates)) (map :gp-mean iterates)))
(ic/view (ip/histogram (map :gp-mean iterates) :nbins 100))
(ic/view (ip/histogram (map :gp-var iterates) :nbins 100))
(ic/view (ip/histogram (map :gp-time-scale iterates) :nbins 100))


write functions generically and then wrap them to modify elements of map
