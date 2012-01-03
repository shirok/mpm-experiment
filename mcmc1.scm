;; -*- coding:utf-8 -*-

;; 1-dimensional Metropolis-Hastings

;; This follows the example shown in Andrieu, C., De Freitas, N., Doucet, A.,
;; and Jordan, M.I., An Introduction to MCMC for Machine Learning,
;; in Machine Learning, 50, 5--43, 2003.

(use gauche.dictionary)
(use srfi-27)
(use srfi-42)

(add-load-path ".")

(use mcmc-util)

;; Target distribution function
(define (p x)
  (+ (* 0.3 (exp (* -0.2 (expt x 2))))
     (* 0.7 (exp (* -0.2 (expt (- x 10) 2))))))

;; Pick a sample x* by proposal distribution - using N(x, 100)
(define (pick1 x) (dpick x 10 -10 20 0.5))

;; Calculate posteriori probability q(x*|x)
(define (q x* x) (/ (φ (/ (- x* x) 10.)) 10.))

;; Single MH step
(define (step x)
  (let ([u  (random-real)]
        [x* (pick1 x)])
    (if (< u (min 1 (/ (* (p x*) (q x* x)) (* (p x) (q x x*)))))
      x*
      x)))

;; Start from x0, run N steps of markov chain
(define (run x0 N)
  (let1 results (make-hash-table 'eqv?)
    (do ([i 0 (+ i 1)]
         [x x0 (step x)])
        [(= i N) results]
      (hash-table-update! results x (cut + 1 <>) 0))))

(define (mcmc1 :key (samples 50000) (file "mcmc1.dat"))
  (with-output-to-file file
    (^[] (let1 result (run 0 samples)
           (do-ec (: x -10 20 0.5)
                  (let1 y (hash-table-get result x 0)
                    (format #t "~a\t~a\n" x (/. y samples))))))))

(define plot
  (^[cmd . args]
    (ecase cmd
      [(close)  (gnuplot 'close)]
      [(result) (let1 p (gnuplot)
                  (format p "set xrange [-10:20]\n")
                  (format p "p(x) = 0.3*exp(-0.2*x*x)+0.7*exp(-0.2*(x-10)*(x-10))\n")
                  (format p "plot p(x), 'mcmc1.dat' using 1:($2*7.5) with boxes\n"))]
      )))


                  
         

    
    


