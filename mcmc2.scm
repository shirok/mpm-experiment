;; -*- coding:utf-8 -*-

;; 2-dimensional Metropolis-Hastings

(use gauche.dictionary)
(use srfi-27)
(use srfi-42)

(add-load-path ".")

(use mcmc-util)

;; Target distribution function
(define (p x y)
  (+ (* 0.3 (exp (+ (* -0.2 (expt (- x 0) 2))
                    (* -0.2 (expt (- y 5) 2)))))
     (* 0.2 (exp (+ (* -0.2 (expt (- x 15) 2))
                    (* -0.2 (expt (- y 10) 2)))))
     (* 0.5 (exp (+ (* -0.1 (expt (- x 10) 2))
                    (* -0.1 (expt (- y 0) 2)))))))

;; Pick a sample (x* y*) by proposal distribution
(define (pick1 x y)
  (values (dpick x 10 -10 20 1)
          (dpick y 10 -10 20 1)))

;; Posteriori probability q([x*,y*]|[x,y])
(define (q x* y* x y)
  (* (/ (φ (/ (- x* x) 10.)) 10.)
     (/ (φ (/ (- y* y) 10.)) 10.)))

;; Single MH step
(define (step x y)
  (let1 u (random-real)
    (receive (x* y*) (pick1 x y)
      (if (< u (min 1 (/ (* (p x* y*) (q x* y* x y))
                         (* (p x y) (q x y x* y*)))))
        (values x* y*)
        (values x y)))))

(define (run x0 y0 N)
  (rlet1 results (make-hash-table 'equal?)
    (let loop ([i 0] [x x0] [y y0])
      (unless (= i N)
        (receive (x y) (step x y)
          (hash-table-update! results (cons x y) (cut + 1 <>) 0)
          (loop (+ i 1) x y))))))

(define (mcmc2 :optional (samples 200000) (file "mcmc2.dat"))
  (with-output-to-file file
    (^[] (let1 result (run 0 0 samples)
           (do-ec (: y -10 20 1)
                  (: x -10 20 1)
                  (begin
                    (when (= x -10) (newline))
                    (let1 z (hash-table-get result (cons x y) 0)
                      (format #t "~a\t~a\t~a\n" x y (/. z samples)))))))))

(define plot
  (^[cmd . args]
    (ecase cmd
      [(close)  (gnuplot 'close)]
      [(result) (let1 p (gnuplot)
                  (format p "set hidden3d\n")
                  (format p "splot 'mcmc2.dat' with lines\n"))])))

  
