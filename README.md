# chicken-lapack

Chicken Scheme bindings for the LAPACK and ATLAS libraries.

ATLAS (http://math-atlas.sourceforge.net) stands for Automatically
Tuned Linear Algebra Software. Its purpose is to provide portably
optimal linear algebra routines. The current version provides a
complete BLAS (http://www.netlib.org/blas) API (for both C and
Fortran77), and a small subset of the LAPACK
(http://www.netlib.org/lapack) API. 

## Naming conventions for routines

Every routine in the LAPACK library comes in four flavors, each
prefixed by the letters S, D, C, and Z, respectively. Each letter
indicates the format of input data:

* S stands for single-precision (32-bit IEEE floating point numbers), 
* D stands for double-precision (64-bit IEEE floating point numbers), 
* C stands for complex numbers (represented by pairs of 32-bit IEEE floating point numbers), 
* Z stands for double complex numbers (represented by pairs of 64-bit IEEE floating point numbers)


In addition, each ATLAS-LAPACK routine in this egg comes in three flavors: 
* Safe, pure: safe routines check the sizes of their input arguments. For example, if a routine is supplied arguments that indicate that an input matrix is of dimensions M-by-N, then the argument corresponding to that matrix is checked that it is of size M * N.  ''Pure'' routines do not alter their arguments in any way. A new matrix or vector is allocated for the return value of the routine. 
* Safe, destructive (suffix: !): safe routines check the sizes of their input arguments. For example, if a routine is supplied arguments that indicate that an input matrix is of dimensions M-by-N, then the argument corresponding to that matrix is checked that it is of size M * N. Destructive routines can modify some or all of their arguments. They  are given names ending in exclamation mark. The matrix factorization routines in LAPACK overwrite the input matrix argument with the result of the factorization, and the linear system solvers overwrite the right-hand side vector with the system solution. Please consult the LAPACK documentation to determine which functions modify their input arguments. 
* Unsafe, destructive (prefix: unsafe-:, suffix: !). Unsafe routines do not check the sizes of their input arguments. They invoke the corresponding ATLAS-LAPACK routines directly. Unsafe routines do not have pure variants. 


For example, function xGESV (N-by-N linear system solver) comes in the following variants: 

<table><tr><th>LAPACK name</th><th>Safe, pure</th><th>Safe, destructive</th><th>Unsafe, destructive</th></tr>
<tr><td>SGESV</td><td>sgesv</td><td>sgesv!</td><td>unsafe-sgesv!</td></tr>
<tr><td>DGESV</td><td>dgesv</td><td>dgesv!</td><td>unsafe-dgesv!</td></tr>
<tr><td>CGESV</td><td>cgesv</td><td>cgesv!</td><td>unsafe-cgesv!</td></tr>
<tr><td>ZGESV</td><td>zgesv</td><td>zgesv!</td><td>unsafe-zgesv!</td></tr>
</table>


## LAPACK datatype conventions

The datatype constructors described below are defined in the
[blas](http://wiki.call-cc.org/eggref/4/blas) library.  The `blas`
library must be loaded before using any of the routines in this
library.

Argument `ORDER` is one of `RowMajor` or `ColMajor` to
indicate that the input and output matrices are in row-major or
column-major form, respectively.

Where present, argument `TRANS` can be one of `NoTrans` or `Trans` to
indicate whether the input matrix is to be transposed or not.

Where present, argument `UPLO` can be one of `Upper` or
`Lower` to indicate whether the upper or lower triangular part
of an input symmetric matrix is to referenced,or to specify the type
of an input triangular matrix.

Where present, argument `DIAG` can be one of `NonUnit` or
`Unit` to indicate whether an input triangular matrix is unit
triangular or not.


## LAPACK driver routines


### General linear system solving


* sgesv:: ORDER * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR * F32VECTOR * S32VECTOR
* dgesv:: ORDER * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR * F64VECTOR * S32VECTOR
+ cgesv:: ORDER * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR * F32VECTOR * S32VECTOR
* zgesv:: ORDER * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR * F64VECTOR * S32VECTOR

The routines compute the solution to a system of linear equations A
* X = B, where A is an N-by-N matrix and X and B are
N-by-NRHS matrices. Optional arguments LDA and LDB are the
leading dimensions of arrays A and B, respectively. LU
decomposition with partial pivoting and row interchanges is used to
factor A as A = P * L * U, where P is a permutation
matrix, L is unit lower triangular, and U is upper
triangular. The factored form of A is then used to solve the
system. The return values are:

* a matrix containing the factors L and U from the factorization A = P*L*U; 
* the N-by-NRHS solution matrix X
*  a vector with pivot indices:  for 1 <= i <= min(M,N), row i of the matrix A was interchanged with row pivot(i)


### Symmetric positive definite linear system solving


<procedure>sposv:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR * F32VECTOR</procedure>

<procedure>dposv:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR * F64VECTOR</procedure>

<procedure>cposv:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR * F32VECTOR</procedure>

<procedure>zposv:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR * F64VECTOR</procedure>



The routines compute the solution to a system of linear equations A * X = B, where A is an N-by-N symmetric positive definite matrix
and X and B are N-by-NRHS matrices. Optional arguments LDA and LDB are the leading dimensions of arrays A and B,
respectively. Cholesky decomposition is used to factor A as * A = U**T * U     if UPLO = 'Upper'
* A = L  * L**T     if UPLO = 'Lower' where U is an upper triangular, and L is a lower triangular matrix. 
The factored form of A is then used to solve the system. The return values are: 
# the factor U or Lfrom the Cholesky factorization, depending on the value of argument UPLO.
# the N-by-NRHS solution matrix X


## LAPACK computational routines


### General matrix factorization


<procedure>sgetrf:: ORDER * M * N * A * [LDA] -> F32VECTOR * S32VECTOR</procedure>

<procedure>dgetrf:: ORDER * M * N * A * [LDA] -> F64VECTOR * S32VECTOR</procedure>

<procedure>cgetrf:: ORDER * M * N * A * [LDA] -> F32VECTOR * S32VECTOR</procedure>

<procedure>zgetrf:: ORDER * M * N * A * [LDA] -> F64VECTOR * S32VECTOR</procedure>



These routines compute an LU factorization of a general M-by-N matrix
A using partial pivoting with row interchanges. Optional argument
LDA is the leading dimension of array A. The return values
are:

* a matrix containing the factors L and U from the factorization A = P*L*U; 
* a vector with pivot indices:  for 1 <= i <= min(M,N), row i of the matrix was interchanged with row pivot(i)


### General linear system solving using factorization


<procedure>sgetrs:: ORDER * TRANSPOSE * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR</procedure>

<procedure>dgetrs:: ORDER * TRANSPOSE * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR</procedure>

<procedure>cgetrs:: ORDER * TRANSPOSE * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR</procedure>

<procedure>zgetrs:: ORDER * TRANSPOSE * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR</procedure>



These routines solve a system of linear equations A * X = B or
A' * X = B with a general N-by-N matrix A using the LU
factorization computed by the xGETRF routines. Argument NRHS is
the number of right-hand sides (i.e. number of columns in
B). Optional arguments LDA and LDB are the leading
dimensions of arrays A and B, respectively. The return value
is the solution matrix X.



### General matrix invert using factorization


<procedure>sgetri:: ORDER * N * A * PIVOT * [LDA] -> F32VECTOR</procedure>

<procedure>dgetri:: ORDER * N * A * PIVOT * [LDA] -> F64VECTOR</procedure>

<procedure>cgetri:: ORDER * N * A * PIVOT * [LDA] -> F32VECTOR</procedure>

<procedure>zgetri:: ORDER * N * A * PIVOT * [LDA] -> F64VECTOR</procedure>



These routines compute the inverse of a matrix using the LU
factorization computed by the xGETRF routines. Argument A must
contain the factors L and U from the LU factorization computed by
xGETRF. Argument PIVOT must be the pivot vector returned by the
factorization routine. Optional argument LDA is the leading
dimension of array A. The return value is the inverse of the
original matrix A.



### Symmetric positive definite matrix factorization


<procedure>spotrf:: ORDER * UPLO * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>dpotrf:: ORDER * UPLO * N * A * [LDA] -> F64VECTOR</procedure>

<procedure>cpotrf:: ORDER * UPLO * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>zpotrf:: ORDER * UPLO * N * A * [LDA] -> F64VECTOR</procedure>



These routines compute the Cholesky factorization of a symmetric positive definite matrix A. The factorization has the form: 

* A = U**T * U     if UPLO = 'Upper'
* A = L  * L**T     if UPLO = 'Lower'

where U is an upper triangular, and L is a lower triangular
matrix. Optional argument LDA is the leading dimension of array
A. The return value is the factor U or Lfrom the Cholesky
factorization, depending on the value of argument UPLO.



### Symmetric positive definite matrix solving using factorization


<procedure>spotrs:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR</procedure>

<procedure>dpotrs:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR</procedure>

<procedure>cpotrs:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F32VECTOR</procedure>

<procedure>zpotrs:: ORDER * UPLO * N * NRHS * A * B * [LDA] * [LDB] -> F64VECTOR</procedure>



These routines solve a system of linear equations A * X = B with a
symmetric positive definite matrix A using the Cholesky
factorization computed by the xPOTRF routines. Argument A is the
triangular factor U or L as computed by xPOTRF. Argument
NRHS is the number of right-hand sides (i.e. number of columns in
B). Argument UPLO indicates whether upper or lower triangle of A
is stored ('Upper' or 'Lower'). Optional arguments LDA and
LDB are the leading dimensions of arrays A and B,
respectively. The return value is the solution matrix X.



### Symmetric positive definite matrix invert using factorization


<procedure>spotri:: ORDER * UPLO * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>dpotri:: ORDER * UPLO * N * A * [LDA] -> F64VECTOR</procedure>

<procedure>cpotri:: ORDER * UPLO * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>zpotri:: ORDER * UPLO * N * A * [LDA] -> F64VECTOR</procedure>



These routines compute the inverse of a symmetric positive definite
matrix A using the Cholesky factorization A = U**T*U or A =
L*L**T computed by xPOTRF. Argument A is the triangular factor
U or L as computed by xPOTRF. Argument UPLO indicates whether
upper or lower triangle of A is stored ('Upper' or
'Lower'). Optional argument LDA is the leading dimension of
array A. The return value is the upper or lower triangle of the
inverse of A.



### Triangular matrix invert


<procedure>strtri:: ORDER * UPLO * DIAG * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>dtrtri:: ORDER * UPLO * DIAG * N * A * [LDA] -> F64VECTOR</procedure>

<procedure>ctrtri:: ORDER * UPLO * DIAG * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>ztrtri:: ORDER * UPLO * DIAG * N * A * [LDA] -> F64VECTOR</procedure>



These routines compute the inverse of a triangular matrix
A. Argument A is the triangular factor U or L as
computed by xPOTRF. Argument UPLO indicates whether upper or lower
triangle of A is stored ('Upper' or 'Lower'). Argument DIAG
indicates whether A is non-unit triangular or unit triangular
('NonUnit' or 'Unit'). Optional argument LDA is the
leading dimension of array A. The return value is the triangular
inverse of the input matrix, in the same storage format.



### Auxilliary routines


<procedure>slauum:: ORDER * UPLO * DIAG * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>dlauum:: ORDER * UPLO * DIAG * N * A * [LDA] -> F64VECTOR</procedure>

<procedure>clauum:: ORDER * UPLO * DIAG * N * A * [LDA] -> F32VECTOR</procedure>

<procedure>zlauum:: ORDER * UPLO * DIAG * N * A * [LDA] -> F64VECTOR</procedure>



These routines compute the product U * U' or L' * L, where the
triangular factor U or L is stored in the upper or lower
triangular part of the array A. Argument UPLO indicates whether
upper or lower triangle of A is stored ('Upper' or
Lower'). Optional argument LDA is the leading dimension of
array A. The return value is the lower triangle of the lower
triangular product, or the upper triangle of upper triangular product,
in the respective storage format.



## Examples
 
```scheme
 (use srfi-4 blas atlas-lapack)
 
 (define order ColMajor)
 (define n 4)
 (define nrhs 1)
 
 ;; Solve the equations
 ;;
 ;; Ax = b,
 ;;
 ;; where A is the general matrix
 
 ;; column-major order
 (define A (f64vector 1.8   5.25   1.58 -1.11  
 		     2.88  -2.95 -2.69 -0.66 
 		     2.05  -0.95 -2.90 -0.59 
 		     -0.89 -3.80 -1.04  0.80))
 ;;
 ;; and b is
 ;;
 (define b (f64vector 9.52 24.35 0.77 -6.22))
 
 ;; A and b are not modified
 (define-values (LU x piv) (dgesv order n nrhs A b))
 
 ;; A is overwritten with its LU decomposition, and 
 ;; b is overwritten with the solution of the system
 (dgesv! order n nrhs A b)
``` 

## Version history

3.0 : Using bind instead of easyffi
2.1 : Ensure that unit test script exists with a non-zero exit status on error (thanks to mario)
2.0 : Eliminated reduntant atlas-lapack: prefix from names of exported symbols

## Requirements

* [blas](http://wiki.call-cc.org/eggref/4/blas)

## License

>
> Copyright 2007-2015 Ivan Raikov
> 
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU General Public License as published by
> the Free Software Foundation, either version 3 of the License, or (at
> your option) any later version.
> 
> This program is distributed in the hope that it will be useful, but
> WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
> General Public License for more details.
> 
> A full copy of the GPL license can be found at
> <http://www.gnu.org/licenses/>.
