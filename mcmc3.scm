;; -*- coding:utf-8 -*-

;; Using MCMC for fitting
;;
;; We have a set of sample points in 2D, and try to find a polynomial
;; that fits them best.
;; First, we assume the dimensions of the search space and run MCMC
;; within it (later we'll include selection of dimensions within the
;; search by using RJMCMC).

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

;; Sample points D = {(x,y)} is given in a format ((x0 . y0) (x1 . y1) ...)
;; Within an area 0 <= x <= 10 and 0 <= y <= 10

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

(define-constant σ 0.002)

;; choose kernel
;; Pick a sample θ* in n-dimensional space
(define (pickN θ) (map (^t (cpick t σ)) θ))

;; Calculate posteriori probability q(θ*|θ)
(define (q θ* θ)
  (apply * (map (^[t* t] (/ (φ (/ (- t* t) σ)) σ)) θ* θ)))

;; Single MH step.  Returns next θ and error^2 value.
(define (step D θ)
  (let* ([u  (random-real)]
         [θ* (pickN θ)]
         [θ1 (if (< u (min 1 (/ (* (p D θ*) (q θ* θ)) (* (p D θ) (q θ θ*)))))
               θ*
               θ)]
         [e  (err D (polynomial θ1))])
    (values θ1 e)))

(define (run D θ0 N)
  (receive (θ e) (step D θ0)
    (let loop ([i 0] [θ θ] [θ-best θ] [e-best e])
      (if (= i N)
        (values θ-best e-best)
        (receive (θ1 e1) (step D θ)
          (print i "\t" e1)
          (if (< e1 e-best)
            (loop (+ i 1) θ1 θ1 e1)
            (loop (+ i 1) θ1 θ-best e-best)))))))

(define (mcmc3 θ0 :key (samples 50000)
                       (datafile "mcmc3.input")
                       (errfile "/dev/null"))
  (let1 D (map (cut apply cons <>) (slices (file->sexp-list datafile) 2))
    (with-output-to-file errfile
      (^[] (run D θ0 samples)))))

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
        (format p "plot 'mcmc3.err' with lines\n")))
    (^[cmd . args]
      (ecase cmd
        [(close)  (gnuplot'close)]
        [(result) (apply show-result args)]
        [(error)  (apply show-error args)]))))

;;
;; some stuff to generate sample inputs
;;
(define (gen-input-1)
  (with-output-to-file "mcmc3.input"
    (^[] (dotimes (i 100)
           (let* ([x (* (random-real) 10)]
                  [y (+ 5 (/ (* (- x 5) (- x 4) (- x 6)) 20)
                        (cpick 0 0.5 -5 5))])
             (when (<= 0 y 10)
               (print x " " y)))))))
