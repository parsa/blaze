//=================================================================================================
/*!
//  \file blaze/math/expressions/SparseMatrix.h
//  \brief Header file for the SparseMatrix base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SPARSEMATRIX_H_
#define _BLAZE_MATH_EXPRESSIONS_SPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup sparse_matrix Sparse Matrices
// \ingroup matrix
*/
/*!\defgroup sparse_matrix_expression Expressions
// \ingroup sparse_matrix
*/
/*!\brief Base class for sparse matrices.
// \ingroup sparse_matrix
//
// The SparseMatrix class is a base class for all sparse matrix classes. It provides an
// abstraction from the actual type of the sparse matrix, but enables a conversion back
// to this type via the Matrix base class.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
class SparseMatrix
   : public Matrix<MT,SO>
{
 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   SparseMatrix() = default;
   SparseMatrix( const SparseMatrix& ) = default;
   SparseMatrix( SparseMatrix&& ) = default;
   ~SparseMatrix() = default;
   SparseMatrix& operator=( const SparseMatrix& ) = default;
   SparseMatrix& operator=( SparseMatrix&& ) = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for row-major sparse matrices.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given row-major sparse
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetLower_backend( SparseMatrix<MT,false>& dm ) -> DisableIf_t< IsUpper_v<MT> >
{
   const size_t m( (~dm).rows() );

   for( size_t i=1UL; i<m; ++i ) {
      (~dm).erase( i, (~dm).begin( i ), (~dm).lowerBound( i, i ) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for column-major sparse matrices.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given column-major sparse
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetLower_backend( SparseMatrix<MT,true>& dm ) -> DisableIf_t< IsUpper_v<MT> >
{
   const size_t m   ( (~dm).rows()    );
   const size_t n   ( (~dm).columns() );
   const size_t jend( min( m, n ) );

   for( size_t j=0UL; j<jend; ++j ) {
      (~dm).erase( j, (~dm).lowerBound( j+1UL, j ), (~dm).end( j ) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetLower() function for lower sparse matrices.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given lower sparse matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline auto resetLower_backend( SparseMatrix<MT,SO>& dm ) -> EnableIf_t< IsUpper_v<MT> >
{
   MAYBE_UNUSED( dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the lower part of the given sparse matrix.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the lower part (excluding the diagonal) of the given sparse matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline void resetLower( SparseMatrix<MT,SO>& dm )
{
   resetLower_backend( ~dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for row-major sparse matrices.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given row-major sparse
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetUpper_backend( SparseMatrix<MT,false>& dm ) -> DisableIf_t< IsLower_v<MT> >
{
   const size_t m   ( (~dm).rows()    );
   const size_t n   ( (~dm).columns() );
   const size_t iend( min( m, n ) );

   for( size_t i=0UL; i<iend; ++i ) {
      (~dm).erase( i, (~dm).lowerBound( i, i+1UL ), (~dm).end( i ) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for column-major sparse matrices.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given column-major sparse
// matrix.
*/
template< typename MT >  // Type of the matrix
inline auto resetUpper_backend( SparseMatrix<MT,true>& dm ) -> DisableIf_t< IsLower_v<MT> >
{
   const size_t n( (~dm).columns() );

   for( size_t j=1UL; j<n; ++j ) {
      (~dm).erase( j, (~dm).begin( j ), (~dm).lowerBound( j, j ) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resetUpper() function for upper sparse matrices.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given upper sparse matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline auto resetUpper_backend( SparseMatrix<MT,SO>& dm ) -> EnableIf_t< IsLower_v<MT> >
{
   MAYBE_UNUSED( dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the upper part of the given sparse matrix.
// \ingroup sparse_matrix
//
// \param matrix The given sparse matrix.
// \return void
//
// This function resets the upper part (excluding the diagonal) of the given sparse matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline void resetUpper( SparseMatrix<MT,SO>& dm )
{
   resetUpper_backend( ~dm );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
