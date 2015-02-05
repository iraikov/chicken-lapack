
(use srfi-1 srfi-4 blas atlas-lapack)

(define order ColMajor)
(define n 4)
(define nrhs 1)

 
(define A (f64vector 1.8   5.25   1.58 -1.11  ;; column-major order
 		     2.88  -2.95 -2.69 -0.66 
 		     2.05  -0.95 -2.90 -0.59 
 		     -0.89 -3.80 -1.04  0.80))
(define b (f64vector 9.52 24.35 0.77 -6.22))
 
;; A and b are not modified
(define-values (LU x piv) (dgesv order n nrhs A b))

(print "x = "  x)

(define eps 1e-14) ;; Used for floating point "equality" test

(assert (every (lambda (u v) (fp< (fp- u v)  eps))
	       (f64vector->list x)
	       (list  1. -1. 3. -5. )))


