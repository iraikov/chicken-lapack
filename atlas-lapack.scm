;;
;; Chicken Scheme bindings for the LAPACK routines in the ATLAS
;; library.
;;
;; Copyright 2007-2015 Ivan Raikov
;;
;; This program is free software: you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation, either version 3 of the
;; License, or (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;; General Public License for more details.
;;
;; A full copy of the GPL license can be found at
;; <http://www.gnu.org/licenses/>.
;;

(module atlas-lapack

   (
     dgesv
     cgesv
     zgesv
     sposv
     dposv
     cposv
     zposv
     sgetrf
     dgetrf
     cgetrf
     zgetrf
     sgetrs
     dgetrs
     cgetrs
     zgetrs
     sgetri
     dgetri
     cgetri
     zgetri
     spotrf
     dpotrf
     cpotrf
     zpotrf
     spotrs
     dpotrs
     cpotrs
     zpotrs
     spotri
     dpotri
     cpotri
     zpotri
     strtri
     dtrtri
     ctrtri
     ztrtri
     slauum
     dlauum
     clauum
     zlauum
     sgesv!
     dgesv!
     cgesv!
     zgesv!
     sposv!
     dposv!
     cposv!
     zposv!
     sgetrf!
     dgetrf!
     cgetrf!
     zgetrf!
     sgetrs!
     dgetrs!
     cgetrs!
     zgetrs!
     sgetri!
     dgetri!
     cgetri!
     zgetri!
     spotrf!
     dpotrf!
     cpotrf!
     zpotrf!
     spotrs!
     dpotrs!
     cpotrs!
     zpotrs!
     spotri!
     dpotri!
     cpotri!
     zpotri!
     strtri!
     dtrtri!
     ctrtri!
     ztrtri!
     slauum!
     dlauum!
     clauum!
     zlauum!
     unsafe-sgesv!
     unsafe-dgesv!
     unsafe-cgesv!
     unsafe-zgesv!
     unsafe-sposv!
     unsafe-dposv!
     unsafe-cposv!
     unsafe-zposv!
     unsafe-sgetrf!
     unsafe-dgetrf!
     unsafe-cgetrf!
     unsafe-zgetrf!
     unsafe-sgetrs!
     unsafe-dgetrs!
     unsafe-cgetrs!
     unsafe-zgetrs!
     unsafe-sgetri!
     unsafe-dgetri!
     unsafe-cgetri!
     unsafe-zgetri!
     unsafe-spotrf!
     unsafe-dpotrf!
     unsafe-cpotrf!
     unsafe-zpotrf!
     unsafe-spotrs!
     unsafe-dpotrs!
     unsafe-cpotrs!
     unsafe-zpotrs!
     unsafe-spotri!
     unsafe-dpotri!
     unsafe-cpotri!
     unsafe-zpotri!
     unsafe-strtri!
     unsafe-dtrtri!
     unsafe-ctrtri!
     unsafe-ztrtri!
     unsafe-slauum!
     unsafe-dlauum!
     unsafe-clauum!
     unsafe-zlauum!)

    (import scheme chicken data-structures foreign)
    (require-extension srfi-4 blas bind)


(bind* #<<EOF
typedef float  CCOMPLEX;
typedef double ZCOMPLEX;


/*
 * Enumerated and derived types from cblas.h
 */

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

/*
 * ===========================================================================
 * Prototypes for LAPACK driver routines
 * ===========================================================================
 */

/* Driver routines for linear equations */

int clapack_sgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  float *A, const int lda, int *ipiv,
                  float *B, const int ldb);

int clapack_dgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  double *A, const int lda, int *ipiv,
                  double *B, const int ldb);

int clapack_cgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  const CCOMPLEX *A, const int lda, int *ipiv,
                  const CCOMPLEX *B, const int ldb);

int clapack_zgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  const ZCOMPLEX *A, const int lda, int *ipiv,
                  const ZCOMPLEX *B, const int ldb);

int clapack_sposv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const int N, const int NRHS, float *A, const int lda,
                  float *B, const int ldb);

int clapack_dposv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const int N, const int NRHS, double *A, const int lda,
                  double *B, const int ldb);

int clapack_cposv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const int N, const int NRHS, const CCOMPLEX *A, const int lda,
                  const CCOMPLEX *B, const int ldb);

int clapack_zposv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const int N, const int NRHS, const ZCOMPLEX *A, const int lda,
                  const ZCOMPLEX *B, const int ldb);

/*
 * ===========================================================================
 * Prototypes for LAPACK computational routines
 * ===========================================================================
 */

/* Computational routines for  linear equations */

/* General matrix factorize */

int clapack_sgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   float *A, const int lda, int *ipiv);

int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   double *A, const int lda, int *ipiv);

int clapack_cgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   const CCOMPLEX *A, const int lda, int *ipiv);

int clapack_zgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   const ZCOMPLEX *A, const int lda, int *ipiv);

/* General matrix solve using factorization */

int clapack_sgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
       const int N, const int NRHS, const float *A, const int lda,
       const int *ipiv, float *B, const int ldb);

int clapack_dgetrs
   (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const double *A, const int lda,
    const int *ipiv, double *B, const int ldb);

int clapack_cgetrs
   (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const const CCOMPLEX *A, const int lda,
    const int *ipiv, const CCOMPLEX *B, const int ldb);

int clapack_zgetrs
   (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const const ZCOMPLEX *A, const int lda,
    const int *ipiv, const ZCOMPLEX *B, const int ldb);


/* General matrix invert using factorization */

int clapack_sgetri(const enum CBLAS_ORDER Order, const int N, float *A,
                   const int lda, const int *ipiv);

int clapack_dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
                   const int lda, const int *ipiv);

int clapack_cgetri(const enum CBLAS_ORDER Order, const int N, const CCOMPLEX *A,
                   const int lda, const int *ipiv);

int clapack_zgetri(const enum CBLAS_ORDER Order, const int N, const ZCOMPLEX *A,
                   const int lda, const int *ipiv);

/* Symmetric/hermitian positive definite matrix factorize */

int clapack_spotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, float *A, const int lda);

int clapack_dpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, double *A, const int lda);

int clapack_cpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const CCOMPLEX *A, const int lda);

int clapack_zpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const ZCOMPLEX *A, const int lda);

/* Symmetric/hermitian positive definite matrix solve using factorization */

int clapack_spotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const float *A, const int lda,
                   float *B, const int ldb);

int clapack_dpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const double *A, const int lda,
                   double *B, const int ldb);

int clapack_cpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const const CCOMPLEX *A, const int lda,
                   const CCOMPLEX *B, const int ldb);

int clapack_zpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const const ZCOMPLEX *A, const int lda,
                   const ZCOMPLEX *B, const int ldb);


/* Symmetric/hermitian positive definite matrix invert using factorization */

int clapack_spotri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, float *A, const int lda);

int clapack_dpotri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, double *A, const int lda);

int clapack_cpotri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const CCOMPLEX *A, const int lda);

int clapack_zpotri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const ZCOMPLEX *A, const int lda);


/* Triangular matrix invert */

int clapack_strtri(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_DIAG Diag,const int N, float *A, const int lda);

int clapack_dtrtri(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_DIAG Diag,const int N, double *A, const int lda);

int clapack_ctrtri(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_DIAG Diag,const int N, const CCOMPLEX *A, const int lda);

int clapack_ztrtri(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_DIAG Diag,const int N, const ZCOMPLEX *A, const int lda);

/* Auxilliary routines  */

int clapack_slauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, float *A, const int lda);

int clapack_clauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const CCOMPLEX *A, const int lda);

int clapack_dlauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, double *A, const int lda);

int clapack_zlauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const ZCOMPLEX *A, const int lda);

EOF
)

(define (atlas-lapack:error x . rest)
  (let ((port (open-output-string)))
    (let loop ((objs (if (symbol? x) rest (cons x rest))))
      (if (null? objs)
    (begin
      (newline port)
      (error (if (symbol? x) x 'atlas-lapack)
       (get-output-string port)))
    (begin (display (car objs) port)
     (display " " port)
     (loop (cdr objs)))))))


(define-syntax lapack-wrap
      (lambda (x r c)
  (let* ((fn      (cadr x))
         (cfname  (string->symbol (conc "clapack_" (symbol->string (car fn)))))
         (ret     (caddr x))
         (errs    (cadddr x))
         (vsize   (car (cddddr x)))
         (copy    (cadr (cddddr x)))
         (fname   (string->symbol (conc (if vsize "" "unsafe-")
                (symbol->string (car fn))
                (if copy "" "!"))))
         (args    (reverse (cdr fn)))

         (fsig     (let loop ((args args) (sig 'rest))
         (if (null? args) (cons fname sig)
             (let ((x (car args)))
         (let ((sig (case x
                ((opiv) sig)
                ((lda)  sig)
                ((ldb)  sig)
                (else   (cons x sig)))))
           (loop (cdr args) sig))))))

         (asize           (r 'asize))
         (bsize           (r 'bsize))

         (%define         (r 'define))
         (%begin          (r 'begin))
         (%let            (r 'let))
         (%cond           (r 'cond))
         (%or             (r 'or))
         (%if             (r 'if))
         (%let-optionals  (r 'let-optionals)))

    `(,%define ,fsig

      (,%let-optionals rest ,(if (memq 'ldb fn)
         `((lda ,(if (memq 'm fn) 'm 'n))
           (ldb ,(if (memq 'm fn) 'm 'n)))
         `((lda ,(if (memq 'm fn) 'm 'n))))

     ,(if vsize
    `(,%begin
       (let ((,asize (,vsize a)))
         ,(if (memq 'm fn)
        `(if (< ,asize (fx* m n))
       (atlas-lapack:error ',fname (conc "matrix A is allocated " ,asize " elements "
                 "but given dimensions are " m " by " n)))
        `(if (< ,asize (fx* n n))
       (atlas-lapack:error ',fname (conc "matrix A is allocated " ,asize " elements "
                 "but given dimensions are " n " by " n)))))
      ,(if (memq 'b fn)
     `(let ((,bsize (,vsize b)))
        ,(if (memq 'nrhs fn)
       `(if (< ,bsize (fx* nrhs n))
            (atlas-lapack:error ',fname (conc "matrix B is allocated " ,bsize " elements "
                "but given dimensions are " n " by " nrhs)))
       `(if (< ,bsize (fx* n 1))
            (atlas-lapack:error ,fname (conc "matrix B is allocated " ,bsize " elements "
                     "but given dimensions are " n " by " 1)))))
     `(,%begin)))
    `(,%begin))

     (let ,(let loop ((fn fn) (bnds '()))
       (if (null? fn) bnds
     (let ((x (car fn)))
       (let ((bnds (case x
         ((opiv)  (cons `(opiv (make-s32vector n)) bnds))
         (else    (if (and copy (memq x ret))
                (cons `(,x (,copy ,x)) bnds)
                bnds)))))
         (loop (cdr fn) bnds)))))

       (let ((info (,cfname . ,(cdr fn))))
     (cond ((= info 0) (values . ,ret))
     ((< info 0) (atlas-lapack:error ',fname (,(car errs) info)))
     ((> info 0) (atlas-lapack:error ',fname (,(cadr errs) info)))))))))
  ))

(define-syntax lapack-wrapx
      (lambda (x r c)
  (let* ((fn     (cadr x))
         (ret    (caddr x))
         (errs   (cadddr x)))
    `(begin
       (lapack-wrap ,(cons (string->symbol (conc "s" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs #f #f)
       (lapack-wrap ,(cons (string->symbol (conc "d" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs #f #f)
       (lapack-wrap ,(cons (string->symbol (conc "c" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs #f #f)
       (lapack-wrap ,(cons (string->symbol (conc "z" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs #f #f)

       (lapack-wrap ,(cons (string->symbol (conc "s" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs f32vector-length #f)
       (lapack-wrap ,(cons (string->symbol (conc "d" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs f64vector-length #f)
       (lapack-wrap ,(cons (string->symbol (conc "c" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs (lambda (v) (fx/ (f32vector-length v) 2)) #f)
       (lapack-wrap ,(cons (string->symbol (conc "z" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs (lambda (v) (fx/ (f64vector-length v) 2)) #f)

       (lapack-wrap ,(cons (string->symbol (conc "s" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs f32vector-length  scopy)
       (lapack-wrap ,(cons (string->symbol (conc "d" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs f64vector-length  dcopy)
       (lapack-wrap ,(cons (string->symbol (conc "c" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs (lambda (v) (fx/ (f32vector-length v) 2)) ccopy)
       (lapack-wrap ,(cons (string->symbol (conc "z" (symbol->string (car fn)))) (cdr fn))
        ,ret ,errs (lambda (v) (fx/ (f64vector-length v) 2)) zcopy))))
      )


  (lapack-wrapx (gesv order n nrhs a lda opiv b ldb)
    (a b opiv)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "upper triangular matrix is singular")))


  (lapack-wrapx (posv order uplo n nrhs a lda b ldb)
    (a b)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) (conc "leading minor of order " i
           " of A is not positive definite"))))

  (lapack-wrapx (getrf order m n a lda opiv)
    (a opiv)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "factor U is singular")))

  (lapack-wrapx (getrs order trans n nrhs a lda ipiv b ldb)
    (b)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "unknown error")))

  (lapack-wrapx (getri order n a lda ipiv)
    (a)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "factor U is singular")))

  (lapack-wrapx (potrf order uplo n a lda)
    (a)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) (conc "leading minor of order " i " is not positive definite"))))

  (lapack-wrapx (potrs order uplo n nrhs a lda b ldb)
    (b)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "unknown error")))

  (lapack-wrapx (potri order uplo n  a lda)
    (a)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) (conc "element " "(" i "," i")" " of factor U or L is zero"))))

  (lapack-wrapx (trtri order uplo diag n a lda)
    (a)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "the triangular matrix is singular")))

  (lapack-wrapx (lauum order uplo n a lda)
    (a)
    ((lambda (i) (conc i "-th argument had an illegal value"))
     (lambda (i) "unknown error")))

)
