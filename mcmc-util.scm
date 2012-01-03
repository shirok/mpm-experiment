;; -*- coding:utf-8 -*-
(define-module mcmc-util
  (use math.const)
  (use gauche.process)
  (use srfi-27)
  (use util.match)
  (use file.util)
  (export φ Φ dpick cpick gnuplot))
(select-module mcmc-util)

(define-constant √2π (sqrt (* pi 2)))

;; Standard normal distribution
(define (φ x) (/ (exp (* -.5 x x)) √2π))

;; Normal cumulative distribution function
;; Abramowitz & Stegun (1964) approximation
(define (Φ x)
  (define (Φ+ x)
    (let* ([t (/ (+ 1 (* 0.2316419 x)))]
           [t2 (* t t)]
           [t3 (* t t2)]
           [t4 (* t t3)]
           [t5 (* t t4)])
      (- 1 (* (φ x)
              (+ (* 0.319381530 t)
                 (* -0.356563782 t2)
                 (* 1.781477937 t3)
                 (* -1.821255978 t4)
                 (* 1.330274429 t5))))))
  (cond [(= x 0) .5]
        [(> x 0) (Φ+ x)]
        [else (- 1 (Φ+ (- x)))]))

;; Pick a discrete sample according to a normal distribution N(μ, σ^2),
;; within a range [min, max] and resolution step
(define (dpick μ σ min max step)
  (let1 u (random-real)
    (do ([x* min (+ x* step)])
        [(or (>= x* max)
             (>= (Φ (/ (- x* μ) σ)) u))
         x*]
      )))

;; Pick a sample from a continuous range, according to a
;; normal distribution N(μ, σ^2).  It is equivalent to solve the following
;; equation for x.
;;   F(x) == Φ((x-μ)/σ) - u = 0  where 0 <= u <= 1
;; We can use Newton-Rhapson.
;;   f(x) == dF(x)/dx == (1/σ)φ((x-μ)/σ)
;; So, given an approximation x_i, we can get a refined appoximation value
;; x_{i+1} as
;;   x_{i+1} = x_i - F(x_i)/f(x_i)
;;
;; If u is very close to 0 or 1, x can be very small or very large.  We clip
;; them by 5σ to avoid instability.

(define (cpick μ σ)
  (let ([min (- μ (* 4 σ))]
        [max (+ μ (* 4 σ))]
        [u (random-real)])
    (let loop ([x0 μ])
      (let1 Fx (- (Φ (/ (- x0 μ) σ)) u)
        (if (< (abs Fx) 10e-6)
          x0
          (let* ([fx (/ (φ (/ (- x0 μ) σ)) σ)]
                 [x1 (- x0 (/ Fx fx))])
            (cond [(<= x1 min) min]
                  [(>= x1 max) max]
                  [else (loop x1)])))))))

;; gnuplot util
(define gnuplot
  (let1 p #f
    (match-lambda*
      [() (process-input
           (or p
               (rlet1 z (run-process '(gnuplot -) :input :pipe :wait #f)
                 (set! p z))))]
      [('close) (when p
                  (close-output-port (process-input p))
                  (process-wait p))])))
