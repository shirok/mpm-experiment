;; -*- coding:utf-8 -*-

;; Apply RJMCMC for generative grammar.
;; Talton et al., Metropolis Procedural Modeling, ACM Trans. Graph. 30(2)

(use gauche.dictionary)
(use gauche.sequence)
(use gauche.process)
(use gauche.record)
(use gauche.threads)
(use gauche.uvector)
(use math.const)
(use util.match)
(use util.list)
(use file.util)
(use srfi-27)
(use srfi-42)
(use srfi-43)
(use gl)
(use gl.glut)
(use gl.simple.viewer)

(add-load-path ".")

(use mcmc-util)

;; Grammar
;;  2-D simple branching model (a modified version of the example in
;;  Talton's paper)
;;
;;  V = {X()}
;;  T = {F(),L(),R(),Z,[,]}
;;  Σ = {l}
;;  ω = {X(1)}
;;  P =
;;    X(l) : l <= M ---(1-l/M)--> F(λ)[L(α,β)X(l+1)][R(α,β)X(l+1)]
;;    X(l) : l <= M ----(l/M)---> F(λ) Z
;;     where λ 〜 N(μ_λ, σ_λ), α 〜 N(μ_α, σ_α), β 〜 N(0, σ_β)
;;
;; Tree representation
;;
;;  <tree> : (Tree <term> ...)
;;  <term> : (X l <term>) | (F (λ)) | (L (α) (β)) | (R (α) (β)) | Z | <tree>
;;
;; The cells (λ), (α) and (β) can be destrictively modified during resampling.

(define M 11)            ;max depth
(define μ_λ 0.5)
(define σ_λ 0.2)
(define μ_α (/ pi 6))
(define σ_α (/ pi 18))
(define μ_β 0)
(define σ_β (/ pi 18))

;; The area we consider as a canvas
(define canvas-x0 -3)
(define canvas-x1  3)
(define canvas-y0  0)
(define canvas-y1  6)

(define-record-type Tree #t #t
  (ω)            ; (X 1 <subtree>)
  (λs)           ; #((λ) ...) descriptive parameter cells
  (αs)           ; #((α) ...) descriptive parameter cells
  (βs)           ; #((β) ...) descriptive parameter cells
  (num-params)   ; sum of the lengths of descriptive parameter vectors
  (Xs)           ; ((X _ _) ...)
  (ws)           ; variable selection weight, see below
  (Σw)           ; sum of variable selection weights, see below
  (Lbuf)         ; likelihood calc buffer; see 'Evaluate likelihood' below
  (p)            ; current likelihood
  (max-tree)     ; A tree that scored the maximum likelihood so far
  (max-p)        ; The maximum likelihood scored so far
  )

;;;
;;; Basic Tree manipulations
;;;

(define (initial-tree)
  (let1 X1 (list 'X 1 #f)
    (rlet1 t (make-Tree X1 '#() '#() '#() 0 `(,X1) '(1) '(1) #f 0.0 '() 0.0)
      ;; calculate initial tree
      (derive! t X1))))

;; Derive a subtree starting from the variable X, and update
;; elements of the tree.
(define (derive! tree X)
  (define (derive-1 Xs)
    (if (null? Xs) (finalize!) (grow! (car Xs) (cdr Xs))))
  (define (finalize!)
    (let ([λs '()] [αs '()] [βs '()] [Xs '()])
      (define (rec tree)
        (match tree
          [('F <λ>) (push! λs <λ>)]
          [('X _ subtree) (push! Xs tree) (rec subtree)]
          [((or 'L 'R) <α> <β>) (push! αs <α>) (push! βs <β>)]
          [('Tree nodes ...) (for-each rec nodes)]
          [_ #f]))
      (rec (Tree-ω tree))
      (Tree-λs-set! tree (list->vector λs))
      (Tree-αs-set! tree (list->vector αs))
      (Tree-βs-set! tree (list->vector βs))
      (Tree-num-params-set! tree (+ (vector-length (Tree-λs tree))
                                    (* (vector-length (Tree-αs tree)) 2)))
      (Tree-Xs-set! tree Xs)
      (recalculate-weights! tree)))
  (define (grow! X Xs)
    (let* ([l (cadr X)] [u (random-real)])
      (if (< u (/ l M))
        (let1 <λ.> (list (cpick μ_λ σ_λ))
          (set! (caddr X) `(Tree (F ,<λ.>) Z))
          (derive-1 Xs))
        (let ([<λ.> (list (cpick μ_λ σ_λ))]
              [<α.> (list (cpick μ_α σ_α))]
              [<β.> (list (cpick μ_β σ_β))]
              [XL (list 'X (+ l 1) #f)]
              [XR (list 'X (+ l 1) #f)])
          (set! (caddr X) `(Tree (F ,<λ.>)
                                 (Tree (L ,<α.> ,<β.>) ,XL)
                                 (Tree (R ,<α.> ,<β.>) ,XR)))
          (derive-1 (list* XL XR Xs))))))
  (derive-1 `(,X)))
  
;; ws is a list of selection weight w_i such that variable X_i is selected
;; to rederiver subtree with the probability qτ(X_i) = w_i/Σw.
;; We set w_i so that qτ(X_i) ∝ b(depth(τ) - depth(X_i)) 
(define (recalculate-weights! tree)
  (let1 maxdepth (apply max (map cadr (Tree-Xs tree)))
    (receive (ws sum)
        (map-accum (^[Xi sum]
                     (let1 wi (expt 2 (- maxdepth (cadr Xi)))
                       (values wi (+ wi sum))))
                   0 (Tree-Xs tree))
      (set! (Tree-ws tree) ws)
      (set! (Tree-Σw tree) sum))))

;; Walk down the tree and calls the handler at every twig and leaf
;; with its coordinates.
(define (walk-node node twig leaf)
  (define (rec node level x y θ)
    (match node
      [('F (λ))
       (let ([x. (+ x (* λ (sin θ)))]
             [y. (+ y (* λ (cos θ)))])
         (when twig (twig x y x. y. level))
         (values x. y. θ))]
      ['Z (when leaf (leaf x y level)) (values x y θ)]
      [('X level subtree) (rec subtree level x y θ)]
      [('Tree nodes ...)
       (let loop ([ns nodes] [x. x] [y. y] [θ. θ])
         (match ns
           [() (values x y θ)]
           [(n . ns) (receive (x. y. θ.) (rec n level x. y. θ.)
                       (loop ns x. y. θ.))]))]
      [((and (or 'L 'R) z) (α) (β))
       (values x y (+ θ (if (eq? z 'L) (+ β α) (- β α))))]))
  (rec node 0 0 0 0))

;; Copy the tree structure to save.
(define (save-tree node)
  (define (rec node)
    (match node
      [('F (λ)) `(F (,λ))]
      ['Z 'Z]
      [('X level subtree) `(X ,level ,(rec subtree))]
      [('Tree nodes ...) `(Tree ,@(map rec nodes))]
      [((and (or 'L 'R) z) (α) (β)) `(,z (,α) (,β))]))
  (rec node))

;; Evaluate likelihood.
;; Desirable input is given in a 32x32 image I.  We calculate
;; exp(-∫ (I(u) - clamp(Σ(Zi(u)),0,1))^2 du) as the likelihood, where
;; I(u) is the pixel value of image I at u, Zi(u) is the value
;; of exp(-((z_i - u) / σ_z)) where z_i is the position of the
;; i-th leaf. 
;; Buffer is the working buffer - caller provides it to avoid allocations.

;; Gaussian kernel.
(define σ_z 3)

(define I_size 32)

(define kernel-center (* 2 (round->exact σ_z))) ;up to 2σ
(define kernel-size (- (* 2 kernel-center) 1))
(define I0 (make-f32vector (* I_size I_size) 1.0))

(define kernel
  (rlet1 v (make-f32vector (* kernel-size kernel-size))
    (do-ec (: y kernel-size)
           (: x kernel-size)
           (f32vector-set! v (+ x (* y kernel-size))
                           (exp (- (/ (+ (expt (- y kernel-center -1) 2)
                                         (expt (- x kernel-center -1) 2))
                                      σ_z)))))))

(define (show-kernel kernel)
  (let1 size (sqrt (uvector-length kernel))
    (do-ec (: y size)
           (: x size)
           (begin
             (format #t " ~4d"
                     (round->exact
                      (* (f32vector-ref kernel (+ x (* y size))) 1000)))
             (when (= x (- size 1)) (newline))))))

(define (likelihood-buf tree)
  (or (Tree-Lbuf tree)
      (rlet1 v (make-f32vector (* I_size I_size))
        (set! (Tree-Lbuf tree) v))))

(define (likelihood tree I)             ;I :: <f32vector>
  (define buf (likelihood-buf tree))
  (define penalty 0)  ; if we're outside of canvas, we penalize it a lot.
  (define (roundI a min max)
    (round->exact (clamp (* (- I_size 1) (/. (- a min) (- max min)))
                         0 (- I_size 1))))
  (define (add-Z x y)
    (if (or (< x canvas-x0) (< canvas-x1 x)
            (< y canvas-y0) (< canvas-y1 y))
      (inc! penalty 5.0)
      (let ([ix (roundI x canvas-x0 canvas-x1)]
            [iy (roundI y canvas-y0 canvas-y1)])
        (do-ec (: y kernel-size)
               (: x kernel-size)
               (let ([jy (+ iy y (- kernel-center))]
                     [jx (+ ix x (- kernel-center))])
                 (when (and (< -1 jy I_size)
                            (< -1 jx I_size))
                   (f32vector-set! buf (+ jx (* jy I_size))
                                   (+ (f32vector-ref kernel
                                                     (+ x (* y kernel-size)))
                                      (f32vector-ref buf
                                                     (+ jx (* jy I_size))))))))
        )))
  
  (f32vector-fill! buf 0)
  (walk-node (Tree-ω tree) #f (^[x y level] (add-Z x y)))
  (f32vector-clamp! buf 0 1)
  (f32vector-sub! buf I)
  ;; Just to calculate likelihood, we can do (exp (- (f32vector-dot buf buf)))
  ;; but we want the squared image for visualization.
  (f32vector-mul! buf buf)
  (exp (- (+ penalty (f32vector-dot buf I0)))))

;;;
;;; RJMCMC
;;;

;; For diffusion move, we resample one random descriptive parameter.  It is
;; independent from the previous state, so q(θ*|θ) and q(θ|θ*) are the same.
;; Unlike mcmc1-mcmc4, we destructively modify the tree. 
(define (diffuse! tree I)
  (let* ([p (Tree-p tree)]
         [rollback! (resample! tree)]
         [p* (likelihood tree I)]
         [u (random-real)])
    (if (< u (min 1 (/ p* p)))
      (Tree-p-set! tree p*)             ;accepted
      (rollback!))))                    ;rejected

;; resample one random parameter, modifies the tree, and returns a
;; procedure that rolls back the tree
(define (resample! tree)    
  (let1 k (random-integer (Tree-num-params tree))
    (cond [(< k (vector-length (Tree-λs tree)))
           (let* ([<λ> (vector-ref (Tree-λs tree) k)]
                  [orig-λ (car <λ>)])
             (set-car! <λ> (cpick μ_λ σ_λ))
             (^[] (set-car! <λ> orig-λ)))]
          [(< k (+ (vector-length (Tree-λs tree))
                   (vector-length (Tree-αs tree))))
           (let* ([k (- k (vector-length (Tree-λs tree)))]
                  [<α> (vector-ref (Tree-αs tree) k)]
                  [orig-α (car <α>)])
             (set-car! <α> (cpick μ_α σ_α))
             (^[] (set-car! <α> orig-α)))]
          [else
           (let* ([k (- k
                        (vector-length (Tree-λs tree))
                        (vector-length (Tree-αs tree)))]
                  [<β> (vector-ref (Tree-βs tree) k)]
                  [orig-β (car <β>)])
             (set-car! <β> (cpick μ_β σ_β))
             (^[] (set-car! <β> orig-β)))])))

;; Jump move.  Pick one random variable and rederiver a subtree below it.
;; Again, we destructively modify the tree.
(define (jump! tree I)
  (receive (X w) (pick-variable tree)
    (let ([rollback! (tree-restorer tree X)]
          [p (Tree-p tree)]
          [qτ (/. w (Tree-Σw tree))]
          [Πφτ (subtree-probability tree X)])
      (derive! tree X)
      (let ([p* (likelihood tree I)]
            [qτ* (/. (find-weight tree X) (Tree-Σw tree))]
            [Πφτ* (subtree-probability tree X)])
        (let ([u (random-real)]
              [A (min 1 (/ (* qτ* p* Πφτ*) (* qτ p Πφτ)))])
          (if (< u A)
            (Tree-p-set! tree p*)       ;accepted
            (rollback!)))))))

;; Variable selection.  We pick a subtree from v with a probability
;; qτ(v) ∝ expt(2, depth(τ)-depth(v)).  Returns the picked node
;; and weight.
(define (pick-variable tree)
  (let1 k (random-integer (Tree-Σw tree))
    (let loop ([ws (Tree-ws tree)] [Xs (Tree-Xs tree)] [sum 0])
      (if (or (null? (cdr ws))
              (< k (+ (car ws) sum)))
        (values (car Xs) (car ws))
        (loop (cdr ws) (cdr Xs) (+ sum (car ws)))))))

(define (find-weight tree X)
  (let loop ([ws (Tree-ws tree)] [Xs (Tree-Xs tree)])
    (if (eq? X (car Xs))
      (car ws)
      (loop (cdr ws) (cdr Xs)))))

;; Walk subtree and multiplies probabilities of each descriptive
;; parameters.
(define (subtree-probability tree X)
  (define (φ. x μ σ) (/ (φ (/ (- x μ) σ)) σ))
  (define (rec node Π)
    (match node
      [('F (λ)) (* Π (φ. λ μ_λ σ_λ))]
      ['Z Π]
      [('X level subtree) (rec subtree Π)]
      [('Tree nodes ...)
       (let loop ([ns nodes] [Π Π])
         (match ns
           [() Π]
           [(n . ns) (loop ns (rec n Π))]))]
      [((or 'L 'R) (α) (β))
       (* Π (φ. α μ_α σ_α) (φ. β μ_β σ_β))]))
  (rec X 1.0))

;; Returns a closure that rollback the tree after derive! operation.
(define (tree-restorer tree X)
  (let ([subtree (caddr X)]
        [λs (Tree-λs tree)] [αs (Tree-αs tree)] [βs (Tree-βs tree)]
        [Xs (Tree-Xs tree)] [ws (Tree-ws tree)] [Σw (Tree-Σw tree)]
        [p  (Tree-p tree)])
    (^[]
      (set! (caddr X) subtree)
      (Tree-λs-set! tree λs) (Tree-αs-set! tree αs) (Tree-βs-set! tree βs)
      (Tree-num-params-set! tree (+ (vector-length λs)
                                    (* (vector-length αs) 2)))
      (Tree-Xs-set! tree Xs) (Tree-ws-set! tree ws) (Tree-Σw-set! tree Σw)
      (Tree-p-set! tree p))))

(define (step tree I)
  (let1 u (random-real)
    (if (< u 0.95)
      (diffuse! tree I)
      (jump! tree I))
    (when (> (Tree-p tree) (Tree-max-p tree))
      (Tree-max-p-set! tree (Tree-p tree))
      (Tree-max-tree-set! tree (save-tree (Tree-ω tree))))))

;;;
;;; Visualization
;;;

(define *tree* #f)

(define *show-max* #f)

(define *show-likelihood* 'likelihood) ; or target, #f

(define *running* #t)

(define (draw-tree tree I)
  (gl-push-matrix)
  (gl-load-identity)
  (when *running* (step tree I))
  (if *show-max*
    (do-tree (Tree-max-tree tree))
    (do-tree (Tree-ω tree)))
  (gl-pop-matrix)
  (show-likelihood-image tree I)
  (glut-post-redisplay))

(define (do-tree node)
  (walk-node node
             (^[x y x. y. level]
               (gl-color (/. 140 255) (/. 92 255) (/. 16 255) 0)
               (gl-line-width (clamp (- 9 level)))
               (gl-begin* GL_LINES (gl-vertex x y) (gl-vertex x. y.)))
             (^[x y level]
               (gl-push-matrix)
               (gl-translate x y 0.01)
               (do-flower)
               (gl-pop-matrix))))

(define (show-likelihood-image tree I)
  (case *show-likelihood*
    [(likelihood)
     (let1 p (likelihood tree I)
       (gl-draw-pixels I_size I_size GL_LUMINANCE GL_FLOAT (Tree-Lbuf tree)))]
    [(target)
     (gl-draw-pixels I_size I_size GL_LUMINANCE GL_FLOAT I)]))
 
(define *flower-id* #f)
(define *flower-quad* #f)

(define (do-flower)
  (unless *flower-id*
    (set! *flower-id* (gl-gen-lists 1))
    (set! *flower-quad* (make <glu-quadric>))
    (gl-new-list *flower-id* GL_COMPILE)
    (gl-color 1.0 (/. 114 255) (/. 114 255) 0)
    (dotimes [n 5]
      (gl-translate 0 0.08 0)
      (glu-disk *flower-quad* 0.0 0.06 10 1)
      (gl-translate 0 -0.08 0)
      (gl-rotate 72 0 0 1))
    (gl-color 1.0 (/. 230 255) (/. 20 255) 0)
    (gl-translate 0 0 0.01)
    (glu-disk *flower-quad* 0.0 0.03 8 1)
    (gl-end-list))
  (gl-call-list *flower-id*))

(define (reshape w h)
  (let1 ratio (/ h w)
    (gl-viewport 0 0 w h)
    (gl-matrix-mode GL_PROJECTION)
    (gl-load-identity)
    (glu-ortho-2d (- canvas-x0 .5) (+ canvas-x1 .5)
                  (- canvas-y0 .5) (+ canvas-y1 .5))
    (gl-matrix-mode GL_MODELVIEW)
    (gl-load-identity)
    ))

(define *tree* #f)

(define (run I)
  (thread-start!
   (make-thread
    (^[]
      (simple-viewer-display
       (^[] (when *tree* (draw-tree *tree* I))))
      (simple-viewer-reshape reshape)
      (simple-viewer-grid #f)
      (simple-viewer-axis #f)
      (simple-viewer-set-key! #f
                              #\m (^[x y]
                                    (set! *show-max* (not *show-max*))
                                    (if *show-max*
                                      (print "Showing maximum fit tree")
                                      (print "Showing every iteration")))
                              #\r (^[x y]
                                    (set! *tree* (initial-tree))
                                    (print "Starting over"))
                              #\i (^[x y]
                                    (case *show-likelihood*
                                      [(likelihood)
                                       (set! *show-likelihood* #f)
                                       (print "Hide likelihood image")]
                                      [(#f)
                                       (set! *show-likelihood* 'target)
                                       (print "Show target image")]
                                      [(target)
                                       (set! *show-likelihood* 'likelihood)
                                       (print "Show likelihood image")]))
                              #\space (^[x y]
                                        (set! *running* (not *running*))
                                        (if *running*
                                          (print "Started")
                                          (print "Stopped")))
                              )
      (simple-viewer-window 'metropolis-procedural-modeling)
      (gl-clear-color 1.0 1.0 1.0 0)
      (simple-viewer-run)))))

;; Some sample target images
(define I-rectangle1 (rlet1 v (make-f32vector (* I_size I_size) 0.0)
                      (do-ec (: y I_size)
                             (: x I_size)
                             (when (< (abs (- x (/ I_size 2)))
                                      (/ I_size 6))
                               (f32vector-set! v (+ x (* y I_size)) 1.0)))))

(define I-rectangle2 (rlet1 v (make-f32vector (* I_size I_size) 0.0)
                      (do-ec (: y I_size)
                             (: x I_size)
                             (when (< (abs (- y (/ I_size 2)))
                                      (/ I_size 6))
                               (f32vector-set! v (+ x (* y I_size)) 1.0)))))

(define I-triangle1 (rlet1 v (make-f32vector (* I_size I_size) 0.0)
                      (do-ec (: y I_size)
                             (: x I_size)
                             (when (< (abs (- x (/ I_size 2)))
                                      (- (/ I_size 2) (/. y 2)))
                               (f32vector-set! v (+ x (* y I_size)) 1.0)))))

(define I-triangle2 (rlet1 v (make-f32vector (* I_size I_size) 0.0)
                      (do-ec (: y I_size)
                             (: x I_size)
                             (when (< (abs (- x (/ I_size 2)))
                                      (/. y 2))
                               (f32vector-set! v (+ x (* y I_size)) 1.0)))))

(define I-checker1 (rlet1 v (make-f32vector (* I_size I_size) 1.0)
                     (do-ec (: y I_size)
                            (: x I_size)
                            (if (<= y (/ I_size 2))
                              (when (< x (/ I_size 2))
                                (f32vector-set! v (+ x (* y I_size)) 0.0))
                              (when (>= x (/ I_size 2))
                                (f32vector-set! v (+ x (* y I_size)) 0.0))))))

(define I-checker2 (rlet1 v (make-f32vector (* I_size I_size) 1.0)
                     (do-ec (: y I_size)
                            (: x I_size)
                            (if (<= y (/ I_size 2))
                              (when (>= x (/ I_size 2))
                                (f32vector-set! v (+ x (* y I_size)) 0.0))
                              (when (< x (/ I_size 2))
                                (f32vector-set! v (+ x (* y I_size)) 0.0))))))

(random-source-randomize! default-random-source)
(set! *tree* (initial-tree))
(glut-init '())
