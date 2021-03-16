//=================================================================================================
/*!
//  \file blaze/math/lapack/pstrf.h
//  \brief Header file for the LAPACK Pivoting Cholesky decomposition functions (pstrf)
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_LAPACK_PSTRF_H_
#define _BLAZE_MATH_LAPACK_PSTRF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/Contiguous.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/clapack/pstrf.h>
#include <blaze/util/Assert.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

template< typename T >
struct RemoveComplex
{
  using type = T;
};

template< typename T >
struct RemoveComplex< complex<T> >
{
  using type = T;
};

//=================================================================================================
//
//  LAPACK LLH PIVOTING (CHOLESKY) DECOMPOSITION FUNCTIONS (PsTRF)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LLH Pivoting (Cholesky) decomposition functions (pstrf) */
//@{
template< typename MT, bool SO >
blas_int_t pstrf( DenseMatrix<MT,SO>& A, char uplo, blas_int_t* piv,
          typename RemoveComplex<ElementType_t<MT>>::type tol );
//@}
//*************************************************************************************************


template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline blas_int_t pstrf( DenseMatrix<MT,SO>& A, char uplo,
  blas_int_t* piv,
  typename RemoveComplex<ElementType_t<MT>>::type tol
    )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   using ET = ElementType_t<MT>;

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows()    ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   if( n == 0 ) {
      return;
   }

   if( IsRowMajorMatrix_v<MT> ) {
      ( uplo == 'L' )?( uplo = 'U' ):( uplo = 'L' );
   }

   blas_int_t lwork( n*2 );
   const std::unique_ptr<ET[]> work( new ET[lwork] );

   blaze::blas_int_t rank;
   pstrf( uplo, n, (*A).data(), lda, piv, &rank, tol, work.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for Cholesky decomposition" );

   for(size_t i=0; i < n; i++)
     piv[i]--;  // From fortran 1-based to C 0-based indexing

   return rank;
}


} // namespace blaze

#endif
