//=================================================================================================
/*!
//  \file blaze/math/views/Subvector.h
//  \brief Header file for all restructuring subvector functions
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Vector.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsMatVecMultExpr.h>
#include <blaze/math/typetraits/IsTVecMatMultExpr.h>
#include <blaze/math/typetraits/IsVecAbsExpr.h>
#include <blaze/math/typetraits/IsVecEvalExpr.h>
#include <blaze/math/typetraits/IsVecScalarDivExpr.h>
#include <blaze/math/typetraits/IsVecScalarMultExpr.h>
#include <blaze/math/typetraits/IsVecTransExpr.h>
#include <blaze/math/typetraits/IsVecVecAddExpr.h>
#include <blaze/math/typetraits/IsVecVecMultExpr.h>
#include <blaze/math/typetraits/IsVecVecSubExpr.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector addition.
// \ingroup views
//
// \param vector The constant vector/vector addition.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the addition.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector addition.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecVecAddExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~vector).leftOperand() , index, size ) +
          subvector( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector subtraction.
// \ingroup views
//
// \param vector The constant vector/vector subtraction.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the subtraction.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector subtraction.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecVecSubExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~vector).leftOperand() , index, size ) -
          subvector( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector multiplication.
// \ingroup views
//
// \param vector The constant vector/vector multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector multiplication.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecVecMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~vector).leftOperand() , index, size ) *
          subvector( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given matrix/vector multiplication.
// \ingroup views
//
// \param vector The constant matrix/vector multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// matrix/vector multiplication.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsMatVecMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   typename VT::LeftOperand  left ( (~vector).leftOperand()  );
   typename VT::RightOperand right( (~vector).rightOperand() );

   return submatrix( left, index, 0UL, size, left.columns() ) * right;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/matrix multiplication.
// \ingroup views
//
// \param vector The constant vector/matrix multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/matrix multiplication.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsTVecMatMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   typename VT::LeftOperand  left ( (~vector).leftOperand()  );
   typename VT::RightOperand right( (~vector).rightOperand() );

   return left * submatrix( right, 0UL, index, right.rows(), size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar multiplication.
// \ingroup views
//
// \param vector The constant vector/scalar multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar multiplication.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecScalarMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~vector).leftOperand(), index, size ) * (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar division.
// \ingroup views
//
// \param vector The constant vector/scalar division.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar division.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecScalarDivExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~vector).leftOperand(), index, size ) / (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector abs operation.
// \ingroup views
//
// \param vector The constant vector abs operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the abs operation.
//
// This function returns an expression representing the specified subvector of the given vector
// abs operation.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecAbsExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return abs( subvector( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector evaluation operation.
// \ingroup views
//
// \param vector The constant vector evaluation operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the evaluation operation.
//
// This function returns an expression representing the specified subvector of the given vector
// evaluation operation.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecEvalExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return eval( subvector( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector transpose operation.
// \ingroup views
//
// \param vector The constant vector transpose operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the transpose operation.
//
// This function returns an expression representing the specified subvector of the given vector
// transpose operation.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecTransExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return trans( subvector( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
