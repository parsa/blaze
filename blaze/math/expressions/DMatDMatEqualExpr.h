//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatDMatEqualExpr.h
//  \brief Header file for the dense matrix/dense matrix equality comparison expression
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATDMATEQUALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATDMATEQUALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/system/Blocking.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL BINARY RELATIONAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two row-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< bool RF         // Relaxation flag
        , typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline bool equal( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         if( !equal( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two column-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< bool RF         // Relaxation flag
        , typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline bool equal( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t j=0UL; j<A.columns(); ++j ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         if( !equal( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two dense matrices with different storage order.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< bool RF       // Relaxation flag
        , typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order
inline bool equal( const DenseMatrix<MT1,SO>& lhs, const DenseMatrix<MT2,!SO>& rhs )
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   const size_t rows   ( A.rows() );
   const size_t columns( A.columns() );
   const size_t block  ( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows; ii+=block ) {
      const size_t iend( ( rows < ii+block )?( rows ):( ii+block ) );
      for( size_t jj=0UL; jj<columns; jj+=block ) {
         const size_t jend( ( columns < jj+block )?( columns ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               if( !equal( A(i,j), B(i,j) ) ) return false;
            }
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side matrix for the comparison.
// \param rhs The right-hand side matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline bool operator==( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline bool operator!=( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return !equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
