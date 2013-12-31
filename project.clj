(defproject bio "0.1.0-SNAPSHOT"
  :description "Experiments with biological sequences"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"],
		 [org.clojure/math.numeric-tower "0.0.2"]]
  :main ^:skip-aot bio.core
  :target-path "target/%s"
  :javac-options ["-target" "1.6" "-source" "1.6" "-Xlint:-options"]
  :profiles {:uberjar {:aot :all}})
