//=================================================================================================
/*!
//  \file blaze/math/traits/RowExprTrait.h
//  \brief Header file for the RowExprTrait class template
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_TRAITS_ROWEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_ROWEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/views/Forward.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the expression type type of a row operation.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a row operation. Given the
// dense or sparse matrix type \a MT, the nested type \a Type corresponds to the resulting return
// type. In case the given type is neither a dense nor a sparse matrix type, the resulting data
// type \a Type is set to \a INVALID_TYPE.
*/
template< typename MT >  // Type of the matrix operand
struct RowExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct DenseResult**************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T >
   struct DenseResult { typedef DenseRow<T,IsRowMajorMatrix<T>::value>  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct SparseResult*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T >
   struct SparseResult { typedef SparseRow<T,IsRowMajorMatrix<T>::value>  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference< typename RemoveCV<MT>::Type >::Type  Tmp;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsExpression<MT>
                      , RowExprTrait<Tmp>
                      , typename If< IsDenseMatrix<MT>
                                   , DenseResult<MT>
                                   , typename If< IsSparseMatrix<MT>
                                                , SparseResult<MT>
                                                , Failure
                                                >::Type
                                   >::Type
                      >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_STATIC_ASSERT( !IsExpression<MT>::value || IsConst<MT>::value ||
                        IsVolatile<MT>::value || IsReference<MT>::value );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
