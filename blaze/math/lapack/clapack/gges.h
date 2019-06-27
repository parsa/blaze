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
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort, int (*selectg)(float*, float*, float*), int n, 
                  float* A, int lda, float* B, int ldb, int* sdim, float* alphar, float* alphal, float* beta, 
                  float* vsl, int ldvsl, float* vsr, int ldvsr, float* work, int lwork, int* bwork, int* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( int ) );
#endif

   sgges_( &jobvsl, &jobvsr, &sort, selectg, &n, A, &lda, B, &ldb, sdim, alphar, alphal, beta, 
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
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort, int (*selectg)(double*, double*, double*), int n, 
                  double* A, int lda, double* B, int ldb, int* sdim, double* alphar, double* alphal, double* beta, 
                  double* vsl, int ldvsl, double* vsr, int ldvsr, double* work, int lwork, int* bwork, int* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( int ) );
#endif

   dgges_( &jobvsl, &jobvsr, &sort, selectg, &n, A, &lda, B, &ldb, sdim, alphar, alphal, beta, 
      vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, info,
      blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1) );
}
//*************************************************************************************************


} // namespace blaze

#endif
