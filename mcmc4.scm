;; -*- coding:utf-8 -*-

;; Like mcmc3, we try to fit a polynomial curve to a given sample points.
;; This time we don't assume the degree of the curve; we use RJMCMC
;; to infer both the degree and the parameters.

(use gauche.dictionary)
(use gauche.sequence)
(use gauche.process)
(use util.match)
(use util.list)
(use file.util)
(use srfi-27)
(use srfi-42)

(add-load-path ".")

(use mcmc-util)

;; 'Jump' occurs between adjacent degrees of polynomials.
;;  n+1 -> n : we jump down one degree with a certain probability δ
;;             (If n+1 == 5 we use 2δ, for we won't jump up).
;;  n -> n+1 : We jump up one degree with a probability δ (If n == 1
;;             we use 2δ)
;;
;;  χ1,2 ≡ {a,b,u}   χ2,1 ≡ {a',b',c'} : u 〜 N(0, σ_u)
;;      f_{1->2} : a' = u/5, b' = a/5, c' = b/5
;;                   |0  0  .2|
;;        J_{1->2} = |.2 0  0 | = (1/5)^3
;;                   |0  .2 0 |
;;      f_{2->1} : a = 5b', b = 5c', u = 5a'  (u will be discarded)
;;                   |0 5 0|
;;        J_{2->1} = |0 0 5| = 5^3
;;                   |5 0 0|
;;      q(k=2|k=1) = 2δ   q_{1->2}(|1,a,b) = 1
;;      q(k=1|k=2) = δ    q_{2->1}(u|2,a',b',c') = (1/σ_u)φ(u/σ_u)
;;

;;  χ2,3 ≡ {a,b,c,u}   χ3,2 ≡ {a',b',c',d'} : u 〜 N(0, σ_u)
;;      f_{2->3} : a' = u/5, b' = a/5, c' = b/5, d' = c/5
;;                   |0  0  0  .2|
;;        J_{2->3} = |.2 0  0  0 | = (1/5)^4
;;                   |0  .2 0  0 |
;;                   |0  0  .2 0 |
;;      f_{3->2} : a = 5b', b = 5c', c = 5d', u = 5a'
;;        J_{3->2} = 5^4
;;      q(k=3|k=2) = δ   q_{2->3}(U|2,a,b,c) = 1
;;      q(k=2|k=3) = δ   q_{3->2}(u|3,a',b',c',d') = (1/σ_u)φ(u/σ_u)
;;
;;  χ3,4 ≡ {a,b,c,d,u}   χ4,3 ≡ {a',b',c',d',e'} : u 〜 N(0, σ)
;;      f_{3->4} : a' = u/5, b' = a/5, c' = b/5, d' = c/5, e' = d/5
;;        J_{3->4} = (1/5)^5
;;      f_{4->3} : a = 5b', b = 5c', c = 5d', d = 5e', u = 5a'
;;        J_{4->3} = 5^5
;;      q(k=4|k=3) = δ   q_{3->4}(U|3,a,b,c,d) = 1
;;      q(k=3|k=4) = δ   q_{4->3}(u|4,a',b',c',d',e') = (1/σ_u)φ(u/σ_u)
;;
;;  χ4,5 ≡ {a,b,c,d,e,u}   χ5,4 ≡ {a',b',c',d',e',f'} : u 〜 N(0, σ)
;;      f_{4->5} : a' = u/5, b' = a/5, c' = b/5, d' = c/5, e' = d/5, f' = e/5
;;        J_{4->5} = (1/5)^6
;;      f_{5->4} : a = 5b', b = 5c', c = 5d', d = 5e', e = 5f', u = 5a'
;;        J_{5->4} = 5^6
;;      q(k=5|k=4) = δ   q_{4->5}(U|4,a,b,c,d,e) = 1
;;      q(k=4|k=5) = 2δ  q_{5->4}(u|5,a',b',c',d',e',f') = (1/σ_u)φ(u/σ_u)


;; Calculate likelihood p(D|θ).  We use gaussian distribution.

(define (p D θ)
  (/ (exp (- (/ (err D (polynomial θ)) (* 2 4)))) 2)) ;σ = 2

(define (err D f)
  (fold (^[d s] (+ s (expt (- (cdr d) (f (car d))) 2))) 0 D))

(define (polynomial θ)
  (match θ
    [(a b)         (^x (+ (* a x) (* 5 b)))]
    [(a b c)       (^x (+ (* a x x) (* 5 b x) (* 5 5 c)))]
    [(a b c d)     (^x (+ (* a x x x) (* 5 b x x) (* 5 5 c x) (* 5 5 5 d)))]
    [(a b c d e)   (^x (+ (* a x x x x) (* 5 b x x x) (* 5 5 c x x)
                          (* 5 5 5 d x) (* 5 5 5 5 e)))]
    [(a b c d e f) (^x (+ (* a x x x x x) (* 5 b x x x x) (* 5 5 c x x x)
                          (* 5 5 5 d x x) (* 5 5 5 5 e x) (* 5 5 5 5 5 f)))]))

(define-constant σ 0.002)               ; for parameter walk
(define-constant σ_u 0.001)             ; for new parameter
(define-constant δ 0.0002)              ; jump

;; choose kernel
;; Pick a sample θ* in n-dimensional space
(define (pickN θ) (map (^t (cpick t σ)) θ))

;; Calculate posteriori probability q(θ*|θ)
(define (q θ* θ)
  (apply * (map (^[t* t] (/ (φ (/ (- t* t) σ)) σ)) θ* θ)))

;; Single MH step.  Returns next m, θ, and error^2 value.
(define (step D m θ)
  (let1 u0 (random-real)
    (case m
      [(1) (if (< u0 (* 2 δ))
             (jump↑ D m (+ m 1) θ)
             (walk D m θ))]
      [(5) (if (< u0 (* 2 δ))
             (jump↓ D m (- m 1) θ)
             (walk D m θ))]
      [else (cond [(< u0 δ)       (jump↓ D m (- m 1) θ)]
                  [(< u0 (* 2 δ)) (jump↑ D m (+ m 1) θ)]
                  [else (walk D m θ)])])))

(define (walk D m θ)
  (let* ([u  (random-real)]
         [θ* (pickN θ)]
         [A  (min 1 (/ (* (p D θ*) (q θ* θ)) (* (p D θ) (q θ θ*))))]
         [θ. (if (< u A) θ* θ)]
         [e  (err D (polynomial θ.))])
    (values m θ. e)))

;; jump from n to m, n > m
(define (jump↓ D n m θ)
  (let* ([u  (random-real)]
         [θ* (pickN (map (cut * 5 <>) θ))]
         [q_n->m (/ (φ (/ (* 5 (car θ)) σ_u)) σ_u)]
         [A  (min 1 (* (/ (p D θ*) (p D θ))
                       (if (= n 5) 2 1) ; q(n|m)/q(m|n)
                       (/ q_n->m)
                       (expt 5 (+ n 1))))]) ;Jacobian
    (if (< u A)
      ;; Jump is accepted.  We discard θ*'s first element.
      (values m (cdr θ*) (err D (polynomial (cdr θ*))))
      ;; Jump is rejected.
      (values n θ (err D (polynomial θ))))))

;; jump from n to m, n < m
(define (jump↑ D n m θ)
  (let* ([u  (random-real)]
         [θadj (map (^t (/ t 5)) (cons (cpick 0 σ_u) θ))]
         [θ* (pickN θadj)]
         [q_m->n (/ (φ (/ (car θadj) σ_u)) σ_u)]
         [A  (min 1 (* (/ (p D θ*) (p D θadj))
                       (if (= n 1) 2 1) ; q(n|m)/q(m|n)
                       q_m->n
                       (/ (expt 5 (+ m 1)))))])   ;Jacobian
    (if (< u A)
      ;; Jump is accepted.
      (values m θ* (err D (polynomial θ*)))
      ;; Jump is rejected.
      (values n θ (err D (polynomial θ))))))

(define (run D m0 θ0 N)
  (receive (m θ e) (step D m0 θ0)
    (let loop ([i 0] [m m] [θ θ] [m-best m] [θ-best θ] [e-best e])
      (if (= i N)
        (values m-best θ-best e-best)
        (receive (m1 θ1 e1) (step D m θ)
          (print i "\t" e1 "\t" m1)
          (if (< e1 e-best)
            (loop (+ i 1) m1 θ1 m1 θ1 e1)
            (loop (+ i 1) m1 θ1 m-best θ-best e-best)))))))

(define (mcmc4 θ0 :key (samples 50000)
                       (datafile "mcmc3.input")
                       (errfile "/dev/null"))
  (let1 D (map (cut apply cons <>) (slices (file->sexp-list datafile) 2))
    (with-output-to-file errfile
      (^[] (run D (- (length θ0) 1) θ0 samples)))))

;;
;; plotting
;;
(define plot
  (let ()
    (define (show-result . results)
      (let1 p (gnuplot)
        (display "set xrange [0:10]\n" p)
        (display "set yrange [0:10]\n" p)
        (format p "plot 'mcmc3.input'~a\n"
                (string-join (map fn results) "," 'prefix))))
    (define (fn θ)
      (match θ
        [(a b) (format "~a*x+(~a)" a (* 5 b))]
        [(a b c) (format "~a*x*x+(~a)*x+(~a)" a (* 5 b) (* 5 5 c))]
        [(a b c d) (format "~a*x**3+(~a)*x*x+(~a)*x+(~a)"
                           a (* 5 b) (* 5 5 c) (* 5 5 5 d))]
        [(a b c d e) (format "~a*x**4+(~a)*x**3+(~a)*x*x+(~a)*x+(~a)"
                             a (* 5 b) (* 5 5 c) (* 5 5 5 d) (* 5 5 5 5 e))]
        [(a b c d e f) (format "~a*x**5+(~a)*x**4+(~a)*x**3+(~a)*x*x+(~a)*x+(~a)"
                             a (* 5 b) (* 5 5 c) (* 5 5 5 d) (* 5 5 5 5 e)
                             (* 5 5 5 5 5 f))]))
    (define (show-error :key (xmin 1000) (ymax '*))
      (let1 p (gnuplot)
        (format p "set xrange [~a:*]\n" xmin)
        (format p "set yrange [0:~a]\n" ymax)
        (format p "set y2range [0:6]\n")
        (format p "plot 'mcmc4.err' using 1:2 axes x1y1 with lines, '' using 1:3 axes x1y2 with steps\n")))
    (^[cmd . args]
      (ecase cmd
        [(close)  (gnuplot'close)]
        [(result) (apply show-result args)]
        [(error)  (apply show-error args)]
        [(m)      (apply show-m args)]))))

