(ns bio.core)

;; program to split header into a hash
;; (defn parse-fasta-header)

;; returns a hash
;; :header contains a hash of header items
;; :seq contains a seq of sequence elements
;;   would be nice to make this lazy, but need to keep the file open to do that
(defn fasta-seq
  "return a hash with :header and :seq keys.  The :seq key returns a seq of elements"
  [filename]
  (let [lines (clojure.string/split-lines (slurp filename))]
    {:header (first lines),
     :seq (seq (clojure.string/join (rest lines)))}))

;; returns a hash
;; :header contains header string
;; :seq contains sequence as a string
(defn fasta
  "parse FASTA file, return a hash with :header and :seq keys.  :seq is a string"
  [filename]
  (let [lines (clojure.string/split-lines (slurp filename))]
    {:header (first lines),
     :seq (clojure.string/join (subvec lines 1))}))

;; get DNA version of P42345.fasta

;; Write DNA to prot

;; Write prot to DNA (disambiguate multiple mappings?)

;; try frequencies, group-by, partition-by, partition(-all) 

;; action of restriction enzymes
;; RFLP - Restriction Fragment Length Polymorphisms
;; (frequencies (map count (restrict enz seq)))

;; better seq rep than chars?

;; 
;; align and find gaps

;;; function to process parsing state in a hash
(defn parse
  "called repeatedly using iterate on a hash represenitng parsing state"
  [state]
  (let [[tok & rem] (:rem state)]
    (println tok)
    (if (not rem) (assoc state :rem rem :end true)
	(assoc state :rem rem))))

; (take-while #(not (:end %)) (iterate parse th)) or
; (->> th (iterate parse) (take-while #(not (:end %))))

;; Given a protein sequence and a DNA sequence, cintinue as long as the
;; frequencies beiween triplets and proteins match
;; return all the matches
(defn mseq
  "walk protein seq ps and DNA seq ds until frequencies mismatch"
  [ds ps]
  (let [trp (partition 3 (:seq ds)),
	pseq (:seq ps),
	pfreq (frequencies pseq),
	dfreq (frequencies trp)]
    (loop [dcode {}, ;; mapping triplet to prot
	   rd trp,
	   rp pseq]
      (let [dh (first rd),
	    ph (first rp)]
	(if (not= (dfreq dh) (pfreq ph)) 
	  (do 
	    (println "mismatch at" dh ph)
	    dcode)
	  (recur (assoc dcode dh ph) (rest rd) (rest rp)))))))

;; keep going on mismatch?

;; (def prot (fasta-seq "P42345.fasta"))
;; (def pfreq (frequencies (:seq prot)))
;; (def rpfreq (clojure.set/map-invert pfreq))
;; (def dna (fasta-seq "mTOR.fasta"))
;; (def dfreq (frequencies (partition 3 (:seq dna))))

;; (defn potential-mappings
;;      "given a freq map of triplets and a freq map of prots, returna possible mapping"
;;      [pfreq dfreq]
;;      (let [rp (clojure.set/map-invert pfreq)]
;;        (for [[k v] dfreq :when (rp v)] [k (rp v)])))

;; ([(\G \T \T) \W] [(\G \A \C) \C] [(\G \A \G) \W] [(\C \T \C) \C] [(\G \C \T) \M])
;; no matches with table
;; order backwards?

(def d2pmap ;; codon table
     {\T {\T {\T \F, \C \F, \A \L, \G \L},
	  \C {\T \S, \C \S, \A \S, \G \S},
	  \A {\T \Y, \C \Y, \A \*, \G \*},
	  \G {\T \C, \C \C, \A \*, \G \W}}
      \C {\T {\T \L, \C \L, \A \L, \G \L},
	  \C {\T \P, \C \P, \A \P, \G \P},
	  \A {\T \H, \C \H, \A \Q, \G \Q},
	  \G {\T \R, \C \R, \A \R, \G \R}}
      \A {\T {\T \I, \C \I, \A \I, \G \M},
	  \C {\T \T, \C \T, \A \T, \G \T},
	  \A {\T \N, \C \N, \A \K, \G \K},
	  \G {\T \S, \C \S, \A \R, \G \R}}
      \G {\T {\T \V, \C \V, \A \V, \G \V},
	  \C {\T \A, \C \A, \A \A, \G \A},
	  \A {\T \D, \C \D, \A \E, \G \E},
	  \G {\T \G, \C \G, \A \G, \G \G}}})

(defn d2p
  "given a DNA char triplet, look up the protein char"
  [[one two three]]
  (((d2pmap one) two) three))

;; (def mtor-d2p (map d2p (partition 3 (:seq dna))))

;; map a DNA sequence to a protein sequence
;; string representation
(defn dna2p
  "given a string DNA sequence,  map to to a string protein sequence"
  [ds]
  (apply str (map d2p (partition 3 ds))))

;; switch rep of :seq to vector, use subvec to do alignments
;; need to get to char vectors
;; or use strings and subs (better than a vector of strings)
;; destructuring works with strings, so (d2p "CAG") works !

;; (def ds (fasta "mTOR.fasta"))
;; (def ps (fasta "p42345.fasta"))
;; (def pstr (:seq ps))
;; (def dpstr (dna2p (:seq ds)))

;; reverse complement - reverse seq, swap G and T, A and C
(defn revcmp
  "create reverse complement of sq, reverse seq, swap G and T, A and C"
  [sq]
  (loop [res '()
	 [head & rem] sq]
    (cond 
     (= head \A) (recur (conj res \C) rem)
     (= head \C) (recur (conj res \A) rem)
     (= head \G) (recur (conj res \T) rem)
     (= head \T) (recur (conj res \G) rem)
     :else (apply str res))))

;; test seqs for alignment
(def tseq "ATCGTATCAGAGCTTTCGTCTGATGCTCGATTGTAA")
(def tal1 "GTATCAGAGCTTTCGTCTGATGCTCGATTGTAA")
(def tal2 "ATCGTATCTTTCGTCTGATGCTCGATTGTAA")
(def tal3 "ATCGTAGAGCTTTCGTCTGATGCTCGATTGTAA")
(def tal4 "ATCGTATCAGAGCGTCTGATGCATTGTAA")
(def tal5 "ATCGTAGAGCCGTCTGCTCGATTGTAA")
(def tal6 "GCTTTCGTCTGATGCTCGATTGTAACAGCAGCAGCAG")

;; align match structure
;; ([pos str] [pos str]...)
;; need src position too
(defrecord Alignment [src-pos    ; position in source sequence
		      targ-pos   ; position in target sequence
		      sq])       ; aligned sequence

(defn walign
  "given a word, find the position it aligns in sequence string sq, -1 if no match"
  [word sq]
  (let [pos (.indexOf sq word)]
    (if (> pos -1) [pos word]
	nil)))

(defn targ-align
  "given an Alignment, return a new alignment with it's targ-pos set or nil if no point is found"
  [algn sq]
  (let [pos (.indexOf sq (:sq algn))]
    (if (> pos -1) (assoc algn :targ-pos pos)
	nil)))

(defn expand-align
  "given a match and a seq, attempt to lengthen the match until it no longer matches"
  [[pos partial] full target]
  (let [maxp (count full),
	maxt (count target)]
   (loop [len (inc (count partial)),
	  npart partial]
     (let [tend (+ pos len)]
       (if (and (<= len maxp) (<= tend maxt))
	 (let [new (subs full 0 len)]
	   (if (=  new (subs target pos tend))
	     (recur (inc len) new)
	     [pos npart]))
	 [pos npart])))))

;; expand short word matches
;; combine expanded matches
;; display alignment

;; merge one pair of alignments
;; return nil if not merge-able
(defn merge-align
  "merge a pair of alignments"
  [[pos1 str1] [pos2 str2]]
  (let [max (count str1)]
    (if (and (<= pos2 (+ pos1 max)) (>= pos2 pos1))
      [pos1 (str (subs str1 0 (- pos2 pos1)) str2)]
      nil)))

;; merge overlapping matches
(defn merge-aligns
  "given a list of alignment matches, merge the overlapping ones"
  [aligns]
  (loop [head (first aligns)
	 rem (rest aligns)
	 res []]
    (if (seq rem)
      (let [candidate (first rem),
	    new-rem (rest rem)]
	(if-let [merged-head (merge-align head candidate)]
	  (recur merged-head new-rem res)
	  (recur candidate new-rem (conj res head))))
      (conj res head))))

;; (keep-indexed (fn [i sq] (Alignment. i 0 (apply str sq))) (partition 5 1 tal2))
;; drop matches inside or before existing matches 
;; march matching short words
(defn align-words
  "align successive words of length word-length from candidate against target"
  [target candidate word-length]
  (let [words (map #(apply str %) (partition word-length 1 candidate))]
    (keep #(walign % target) words)))

(defn align
  "create an alignment of candidate against target based on a given word length"
  [target candidate word-length]
  (merge-aligns (align-words target candidate word-length)))

;; bug in merge-aligns on length under 4?
;; short seqs likely to match early in the target, match list will be out of order
;; need to sort on position, eliminate shorter matches, or start newer matches after previous matches
;; for now, just return nil in merge-align if second pos if before first
;; need to keep track of position in candidate as well as target

(defn display-align
  "print alignment given target and alignment result"
  [target aligns]
  )
