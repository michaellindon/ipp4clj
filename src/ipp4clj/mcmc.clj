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
(def gp-var 2)
(def gp-time-scale 0.2)
(def F (sample-gp {} gp-var gp-time-scale))
(def f (comp first F))
(def max-intensity 1000)
(def gp-mean 1)
(def intensity (comp (partial * max-intensity) probit (partial + gp-mean ) f))
(ic/view (ip/function-plot f 0 1))
(def obs (sample-ppp intensity max-intensity 0 1))
(def thin-intensity (fn [x] (- max-intensity (intensity x))))
(def aug (sample-ppp thin-intensity max-intensity 0 1))
(def obs-y (map (fn [x] (sample-right-normal (+ gp-mean (f x)) 1)) obs))
(def aug-y (map (fn [x] (sample-left-normal (+ gp-mean (f x)) 1)) aug))
(def obs-map (into (avl/sorted-map) (map vector obs obs-y)))
(def aug-map (into (avl/sorted-map) (map vector aug aug-y)))
(def foo (merge obs-map aug-map))
(ic/view (ip/scatter-plot (keys foo) (vals foo)))

(def initial-state {:F F
                    :gp-var gp-var
                    :gp-time-scale gp-time-scale
                    :gp-mean gp-mean
                    :max-intensity max-intensity
                    :obs-times obs
                    :aug-times aug
                    :y (merge obs-map aug-map)
                    :start-t 0
                    :end-t 1})

(defn ulogpdf-gamma [a b x]
 (if (neg? x)
     Double/NEGATIVE_INFINITY
     (- (* (dec a) (log x)) (* b x))))

(defn update-gp-time-scale [state]
 (let [{y :y
        gp-var :gp-var
        gp-time-scale :gp-time-scale} state
        times (keys y)
        obs (vals y)
        conditional (fn [x] (+ (ulogpdf-gamma 2 1 x) (log-likelihood times obs 1 gp-var x)))]
   (assoc state :gp-time-scale (sample-slice conditional 1 gp-time-scale))))

(defn update-gp-var [state]
 (let [{y :y
        gp-var :gp-var
        gp-time-scale :gp-time-scale} state
        times (keys y)
        obs (vals y)
        conditional (fn [x] (+ (ulogpdf-gamma 2 1 x) (log-likelihood times obs 1 x gp-time-scale)))]
   (assoc state :gp-var (sample-slice conditional 1 gp-var))))


(defn update-gp-mean [state]
 (let [{y :y
        gp-var :gp-var
        gp-time-scale :gp-time-scale} state
        times (keys y)
        obs (vals y)
        mus2 (mean-conditional times obs 1 gp-var gp-time-scale)
        new-mean (sample-normal 1 :mean (first mus2) :sd (sqrt (second mus2)))]
   (assoc state :gp-mean new-mean)))


(defn update-aug-times [state]
 (let [{F :F
        max-intensity :max-intensity
        start-t :start-t
        end-t :end-t
        gp-mean :gp-mean} state
        f (comp first F)
        intensity (comp (partial * max-intensity) probit (partial + gp-mean) f)
        thin-intensity (fn [x] (- max-intensity (intensity x)))]
   (assoc state :aug-times (sample-ppp thin-intensity max-intensity start-t end-t))))

(defn update-y [state]
 (let [{obs-times :obs-times
        aug-times :aug-times
        gp-mean :gp-mean
        F :F} state
        f (comp first F)
        obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
        aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
        obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
        aug-map (into (avl/sorted-map) (map vector aug-times aug-y))]
   (assoc state :y (merge obs-map aug-map))))

(defn update-F [state]
 (let [{gp-var :gp-var
        gp-time-scale :gp-time-scale
        gp-mean :gp-mean
        y :y} state
        times (keys y)
        obs (map (fn [x] (- x gp-mean)) (vals y))
        pivots (FFBS times obs 1 gp-var gp-time-scale)]
   (assoc state :F (sample-gp pivots gp-var gp-time-scale))))


(def transition (comp
                      update-y
                      update-aug-times
                      update-F
                      update-gp-mean
                      update-gp-var
                      update-gp-time-scale))

(def mcmc (iterate transition initial-state))
(defn mean-function [x] (mean (map (fn [y] (first ((:F y) x))) (take 200 mcmc))))
(defn lower-function [x] (quantile (map (fn [y] (first ((:F y) x))) (take 200 mcmc)) :probs 0.025))
(defn upper-function [x] (quantile (map (fn [y] (first ((:F y) x))) (take 200 mcmc)) :probs 0.975))
(def learned (:F (mget mcmc 100)))
(ic/view (ip/function-plot (comp first learned) 0 1))
(ic/view (ip/function-plot f 0 1))
(ic/view (ip/function-plot mean-function 0 1))
(def learned (:y (transition (transition (transition initial-state)))))




(def ribbon (ip/function-plot mean-function 0 1))
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 0) 0 java.awt.Color/red)
(ip/add-function ribbon upper-function 0 1 :series-label 1)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 1) 0 java.awt.Color/pink)
(ip/add-function ribbon lower-function 0 1 :series-label 2)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 2) 0 java.awt.Color/pink)
(ip/add-function ribbon f 0 1 :series-label 3)
(.setSeriesPaint (.getRenderer (.getPlot ribbon) 3) 0 java.awt.Color/black)
(ic/view ribbon)


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







(def Atrials (take 3 (repeatedly (fn [] {:obs-times (sample-ppp intensity max-intensity 0 1)
                                         :aug-times []
                                         :y (avl/sorted-map)
                                         :start-t 0
                                         :end-t 1}))))



(first (:obs-times (first Atrials)))
(first (:obs-times (second Atrials)))
(type Atrials)
                                        ;Atrials should be a key-val pair in the state map


(defn update-y [state]
  (let [{obs-times :obs-times
         aug-times :aug-times
         gp-mean :gp-mean
         F :F} state
        f (comp first F)
        obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
        aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)
        obs-map (into (avl/sorted-map) (map vector obs-times obs-y))
        aug-map (into (avl/sorted-map) (map vector aug-times aug-y))]
    (assoc state :y (merge obs-map aug-map))))

(defn update-Atrials [state]
  (let [{F :F
         max-intensity :max-intensity
         gp-mean :gp-mean state
         f (comp first F)
         intensity (comp (partial * max-intensity) probit (partial + gp-mean) f)
         thin-intensity (fn [x] (- max-intensity (intensity x)))
         Atrials :Atrials} state
        update-aug-times (fn [trial]
                           (let [{start-t :start-t
                                  end-t :end-t} trial]
                             (assoc :aug-times (sample-ppp thin-intensity max-intensity start-t end-t))))
        update-y-times (fn [trial]
                         (let [{obs-times :obs-times
                                aug-times :aug-times} trial
                               obs-y (map (fn [t] (sample-right-normal (+ gp-mean (f t)) 1 )) obs-times)
                               aug-y (map (fn [t] (sample-left-normal (+ gp-mean (f t)) 1 )) aug-times)]
                           (assoc trial :y (merge obs-map aug-map))))
        update-trials (comp update-y-times update-aug-times)
        ]
    (assoc state :Atrials (map update-trials Atrials))))

;map a function across trials which gives new Atrials vector
