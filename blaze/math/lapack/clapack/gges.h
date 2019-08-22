//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/gees.h
//  \brief Header file for the CLAPACK gges wrapper functions
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_GGES_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_GGES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(INTEL_MKL_VERSION)
extern "C" {

void sgges_( char* jobvsl, char* jobvsr, char* sort, int (*selectg)(float*, float*, float*), int* n, 
   float* A, int* lda, float* B, int* ldb, int* sdim, float* alphar, float* alphal, float* beta, 
   float* vsl, int* ldvsl, float* vsr, int* ldvsr, float* work, int* lwork, int* bwork, int* info,
   blaze::fortran_charlen_t njobvsl, blaze::fortran_charlen_t njobvsr, blaze::fortran_charlen_t nsort );
   
void dgges_( char* jobvsl, char* jobvsr, char* sort, int (*selectg)(double*, double*, double*), int* n, 
   double* A, int* lda, double* B, int* ldb, int* sdim, double* alphar, double* alphal, double* beta, 
   double* vsl, int* ldvsl, double* vsr, int* ldvsr, double* work, int* lwork, int* bwork, int* info,
   blaze::fortran_charlen_t njobvsl, blaze::fortran_charlen_t njobvsr, blaze::fortran_charlen_t nsort );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK GENERALIZED SCHUR DECOMPOSITION FUNCTIONS (GGES)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK generalized Schur decomposition functions (gges) */
//@{
inline void gges( char jobvsl, char jobvsr, char sort, int (*selectg)(float*, float*, float*), int n, 
                  float* A, int lda, float* B, int ldb, int sdim, float* alphar, float* alphal, float* beta, 
                  float* vsl, int ldvsl, float* vsr, int ldvsr, float* work, int lwork, int* bwork, int info );

inline void gges( char jobvsl, char jobvsr, char sort, int (*selectg)(double*, double*, double*), int n, 
                  double* A, int lda, double* B, int ldb, int sdim, double* alphar, double* alphal, double* beta, 
                  double* vsl, int ldvsl, double* vsr, int ldvsr, double* work, int lwork, int* bwork, int info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur decomposition
// for a pair of non-symmetric single precision matrices.
// \ingroup lapack_eigenvalue
//
// GGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),
//  the generalized eigenvalues, the generalized real Schur form (S,T),
//  optionally, the left and/or right matrices of Schur vectors (VSL and
//  VSR). This gives the generalized Schur factorization
// 
//           (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
// 
//  Optionally, it also orders the eigenvalues so that a selected cluster
//  of eigenvalues appears in the leading diagonal blocks of the upper
//  quasi-triangular matrix S and the upper triangular matrix T.The
//  leading columns of VSL and VSR then form an orthonormal basis for the
//  corresponding left and right eigenspaces (deflating subspaces).
//
// For more information on the sgges() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \param jobvsl
//           = 'N':  do not compute the left Schur vectors;
//           = 'V':  compute the left Schur vectors.
// \param jobvsr
//           = 'N':  do not compute the right Schur vectors;
//           = 'V':  compute the right Schur vectors.
// \param sort
//           Specifies whether or not to order the eigenvalues on the
//           diagonal of the generalized Schur form.
//           = 'N':  Eigenvalues are not ordered;
//           = 'S':  Eigenvalues are ordered (see \a selctg);
// \param selctg
//           a pointer to a function of three single precision floating point arguments
//           returning a boolean value.
//           If \a sort = 'N', \a selctg is not referenced.
//           If \a sort = 'S', \a selctg is used to select eigenvalues to sort
//           to the top left of the Schur form.
//           An eigenvalue (alphar[j]+I*alphai[j])/beta[j] is selected if
//           selctg(&alphar[j],&alphai[j],&beta[j]) is true; i.e. if either
//           one of a complex conjugate pair of eigenvalues is selected,
//           then both complex eigenvalues are selected.
// 
//           Note that in the ill-conditioned case, a selected complex
//           eigenvalue may no longer satisfy selctg(&alphar[j],&alphai[j],
//           &beta[j]) = TRUE after ordering. \a info is to be set to \a n+2
//           in this case.
// \param	n	
//           The order of the matrices \a A, \a B, \a VSL, and \a vsr.  \a n >= 0.
// \param	A	
//           On entry, the first of the pair of matrices.
//           On exit, \a A has been overwritten by its generalized Schur
//           form S.
// \param	lda
//           The leading dimension of \a A.  \a lda >= max(1,\a n).
// \param	B	
//           On entry, the second of the pair of matrices.
//           On exit, \a B has been overwritten by its generalized Schur
//           form T.
// \param	ldb	
//           The leading dimension of \a B.  \a ldb >= max(1,\a n).
// \param	sdim	
//           If \a sort = 'N', \a *\a sdim = 0.
//           If \a sort = 'S', *\a sdim = number of eigenvalues (after sorting)
//           for which \a selctg is true.  (Complex conjugate pairs for which
//           \a selctg is true for either eigenvalue count as 2.)
// \param	alphar	
//           real part of eigenvalue numerator, dimension (\a n)
// \param	alphai	
//           imaginary part of eigenvalue numerator, dimension (\a n)
// \param	beta
//           eigenvalue denominator, dimension (\a n).
//           On exit, (alphar[j] + I*alphai[j])/beta[j], j=0,...,n-1, will
//           be the generalized eigenvalues.  alphar[j] + I*alphai[j],
//           and  beta[j],j=0,...,n-1 are the diagonals of the complex Schur
//           form (S,T) that would result if the 2-by-2 diagonal blocks of
//           the real Schur form of (A,B) were further reduced to
//           triangular form using 2-by-2 complex unitary transformations.
//           If alphai[j] is zero, then the j-th eigenvalue is real; if
//           positive, then the j-th and (j+1)-st eigenvalues are a
//           complex conjugate pair, with alphai[j+1] negative.
// 
//           Note: the quotients alphar[j]/beta[j] and alphai[j]/beta[j]
//           may easily over- or underflow, and beta[j] may even be zero.
//           Thus, the user should avoid naively computing the ratio.
//           However, alphar and alphai will be always less than and
//           usually comparable with norm(A) in magnitude, and beta always
//           less than and usually comparable with norm(B).
// \param	vsl	
//           If \a jobvsl = 'V', \a vsl will contain the left Schur vectors.
//           Not referenced if \a jobvsl = 'N'.
// \param	ldvsl
//           The leading dimension of the matrix \a vsl. \a lsvsl >=1, and
//           if \a jobvsl = 'V', \a ldvsl >= \a n.
// \param	vsr
//           If \a jobvsr = 'V', \a vsr will contain the right Schur vectors.
//           Not referenced if \a jobvsr = 'N'.
// \param	ldvsr	
//           The leading dimension of the matrix \a vsr. \a ldvsr >= 1, and
//           if \a jobvsr = 'V', \a ldvsr >= \a n.
// \param	work
//           is a floating point working array, dimension (MAX(1,\a lwork))
//           On exit, if *\a info = 0, \a work[0] returns the optimal \a lwork.
// \param	lwork	
//           The dimension of the array \a work.
//           If \a n = 0, \a lwork >= 1, else \a lwork >= 8*\a n+16.
//           For good performance, \a lwork must generally be larger.
// 
//           If \a lwork = -1, then a workspace query is assumed; the routine
//           only calculates the optimal size of the work array, returns
//           this value as the first entry of the work array, and no error
//           message related to \a lwork is issued by XERBLA.
// \param	bwork	
//           is a logical working array, dimension (\a n)
//           Not referenced if \a sort = 'N'.
// \param	info	
//           = 0:  successful exit
//           < 0:  if \a info = -i, the i-th argument had an illegal value.
//           = 1,...,\a n:
//                 The QZ iteration failed.  (A,B) are not in Schur
//                 form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
//                 be correct for j=\a info,...,\a n-1.
//           > \a n:  =\a n+1: other than QZ iteration failed in DHGEQZ.
//                 =\a n+2: after reordering, roundoff changed values of
//                       some complex eigenvalues so that leading
//                       eigenvalues in the Generalized Schur form no
//                       longer satisfy \a selctg=TRUE.  This could also
//                       be caused due to scaling.
//                 =\a n+3: reordering failed in DTGSEN.
//
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort, int (*selectg)(float*, float*, float*), int n, 
                  float* A, int lda, float* B, int ldb, int* sdim, float* alphar, float* alphai, float* beta, 
                  float* vsl, int ldvsl, float* vsr, int ldvsr, float* work, int lwork, int* bwork, int* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( int ) );
#endif

   sgges_( &jobvsl, &jobvsr, &sort, selectg, &n, A, &lda, B, &ldb, sdim, alphar, alphai, beta, 
      vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, info, 
      blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur decomposition
// for a pair of non-symmetric double precision matrices.
// \ingroup lapack_eigenvalue
//
// GGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),
//  the generalized eigenvalues, the generalized real Schur form (S,T),
//  optionally, the left and/or right matrices of Schur vectors (VSL and
//  VSR). This gives the generalized Schur factorization
// 
//           (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
// 
//  Optionally, it also orders the eigenvalues so that a selected cluster
//  of eigenvalues appears in the leading diagonal blocks of the upper
//  quasi-triangular matrix S and the upper triangular matrix T.The
//  leading columns of VSL and VSR then form an orthonormal basis for the
//  corresponding left and right eigenspaces (deflating subspaces).
//
// For more information on the dgges() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
// 
// \param jobvsl
//           = 'N':  do not compute the left Schur vectors;
//           = 'V':  compute the left Schur vectors.
// \param jobvsr
//           = 'N':  do not compute the right Schur vectors;
//           = 'V':  compute the right Schur vectors.
// \param sort
//           Specifies whether or not to order the eigenvalues on the
//           diagonal of the generalized Schur form.
//           = 'N':  Eigenvalues are not ordered;
//           = 'S':  Eigenvalues are ordered (see \a selctg);
// \param selctg
//           a pointer to a function of three double precision floating point arguments
//           returning a boolean value.
//           If \a sort = 'N', \a selctg is not referenced.
//           If \a sort = 'S', \a selctg is used to select eigenvalues to sort
//           to the top left of the Schur form.
//           An eigenvalue (alphar[j]+I*alphai[j])/beta[j] is selected if
//           selctg(&alphar[j],&alphai[j],&beta[j]) is true; i.e. if either
//           one of a complex conjugate pair of eigenvalues is selected,
//           then both complex eigenvalues are selected.
// 
//           Note that in the ill-conditioned case, a selected complex
//           eigenvalue may no longer satisfy selctg(&alphar[j],&alphai[j],
//           &beta[j]) = TRUE after ordering. \a info is to be set to \a n+2
//           in this case.
// \param	n	
//           The order of the matrices \a A, \a B, \a VSL, and \a vsr.  \a n >= 0.
// \param	A	
//           On entry, the first of the pair of matrices.
//           On exit, \a A has been overwritten by its generalized Schur
//           form S.
// \param	lda
//           The leading dimension of \a A.  \a lda >= max(1,\a n).
// \param	B	
//           On entry, the second of the pair of matrices.
//           On exit, \a B has been overwritten by its generalized Schur
//           form T.
// \param	ldb	
//           The leading dimension of \a B.  \a ldb >= max(1,\a n).
// \param	sdim	
//           If \a sort = 'N', \a *\a sdim = 0.
//           If \a sort = 'S', *\a sdim = number of eigenvalues (after sorting)
//           for which \a selctg is true.  (Complex conjugate pairs for which
//           \a selctg is true for either eigenvalue count as 2.)
// \param	alphar	
//           real part of eigenvalue numerator, dimension (\a n)
// \param	alphai	
//           imaginary part of eigenvalue numerator, dimension (\a n)
// \param	beta
//           eigenvalue denominator, dimension (\a n).
//           On exit, (alphar[j] + I*alphai[j])/beta[j], j=0,...,n-1, will
//           be the generalized eigenvalues.  alphar[j] + I*alphai[j],
//           and  beta[j],j=0,...,n-1 are the diagonals of the complex Schur
//           form (S,T) that would result if the 2-by-2 diagonal blocks of
//           the real Schur form of (A,B) were further reduced to
//           triangular form using 2-by-2 complex unitary transformations.
//           If alphai[j] is zero, then the j-th eigenvalue is real; if
//           positive, then the j-th and (j+1)-st eigenvalues are a
//           complex conjugate pair, with alphai[j+1] negative.
// 
//           Note: the quotients alphar[j]/beta[j] and alphai[j]/beta[j]
//           may easily over- or underflow, and beta[j] may even be zero.
//           Thus, the user should avoid naively computing the ratio.
//           However, alphar and alphai will be always less than and
//           usually comparable with norm(A) in magnitude, and beta always
//           less than and usually comparable with norm(B).
// \param	vsl	
//           If \a jobvsl = 'V', \a vsl will contain the left Schur vectors.
//           Not referenced if \a jobvsl = 'N'.
// \param	ldvsl
//           The leading dimension of the matrix \a vsl. \a lsvsl >=1, and
//           if \a jobvsl = 'V', \a ldvsl >= \a n.
// \param	vsr
//           If \a jobvsr = 'V', \a vsr will contain the right Schur vectors.
//           Not referenced if \a jobvsr = 'N'.
// \param	ldvsr	
//           The leading dimension of the matrix \a vsr. \a ldvsr >= 1, and
//           if \a jobvsr = 'V', \a ldvsr >= \a n.
// \param	work
//           is a floating point working array, dimension (MAX(1,\a lwork))
//           On exit, if *\a info = 0, \a work[0] returns the optimal \a lwork.
// \param	lwork	
//           The dimension of the array \a work.
//           If \a n = 0, \a lwork >= 1, else \a lwork >= 8*\a n+16.
//           For good performance, \a lwork must generally be larger.

//           If \a lwork = -1, then a workspace query is assumed; the routine
//           only calculates the optimal size of the work array, returns
//           this value as the first entry of the work array, and no error
//           message related to \a lwork is issued by XERBLA.
// \param	bwork	
//           is a logical working array, dimension (\a n)
//           Not referenced if \a sort = 'N'.
// \param	info	
//           = 0:  successful exit
//           < 0:  if \a info = -i, the i-th argument had an illegal value.
//           = 1,...,\a n:
//                 The QZ iteration failed.  (A,B) are not in Schur
//                 form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
//                 be correct for j=\a info,...,\a n-1.
//           > \a n:  =\a n+1: other than QZ iteration failed in DHGEQZ.
//                 =\a n+2: after reordering, roundoff changed values of
//                       some complex eigenvalues so that leading
//                       eigenvalues in the Generalized Schur form no
//                       longer satisfy \a selctg=TRUE.  This could also
//                       be caused due to scaling.
//                 =\a n+3: reordering failed in DTGSEN.
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort, int (*selectg)(double*, double*, double*), int n, 
                  double* A, int lda, double* B, int ldb, int* sdim, double* alphar, double* alphai, double* beta, 
                  double* vsl, int ldvsl, double* vsr, int ldvsr, double* work, int lwork, int* bwork, int* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( int ) );
#endif

   dgges_( &jobvsl, &jobvsr, &sort, selectg, &n, A, &lda, B, &ldb, sdim, alphar, alphai, beta, 
      vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, info,
      blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1) );
}
//*************************************************************************************************


} // namespace blaze

#endif
