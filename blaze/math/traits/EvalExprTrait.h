//=================================================================================================
/*!
//  \file blaze/math/traits/EvalExprTrait.h
//  \brief Header file for the EvalExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_EVALEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_EVALEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DMatEvalExprTrait.h>
#include <blaze/math/traits/DVecEvalExprTrait.h>
#include <blaze/math/traits/SMatEvalExprTrait.h>
#include <blaze/math/traits/SVecEvalExprTrait.h>
#include <blaze/math/traits/TDMatEvalExprTrait.h>
#include <blaze/math/traits/TDVecEvalExprTrait.h>
#include <blaze/math/traits/TSMatEvalExprTrait.h>
#include <blaze/math/traits/TSVecEvalExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsNumeric.h>
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
/*!\brief Evaluation of the return type of an evaluation expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of an evaluation expression.
// Given the type \a T, which must be either a vector or matrix type, the nested type \a Type
// corresponds to the resulting return type. In case the type of \a T doesn't fit or if no
// evaluation operation exists for the type, the resulting data type \a Type is set to
// \a INVALID_TYPE.
*/
template< typename T >  // Type of the evaluation operand
struct EvalExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<T>::value || IsVolatile<T>::value || IsReference<T>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsMatrix<T>
                      , typename If< IsDenseMatrix<T>
                                   , typename If< IsRowMajorMatrix<T>
                                                , DMatEvalExprTrait<T>
                                                , TDMatEvalExprTrait<T>
                                                >::Type
                                   , typename If< IsRowMajorMatrix<T>
                                                , SMatEvalExprTrait<T>
                                                , TSMatEvalExprTrait<T>
                                                >::Type
                                   >::Type
                      , typename If< IsVector<T>
                                   , typename If< IsDenseVector<T>
                                                , typename If< IsTransposeVector<T>
                                                             , TDVecEvalExprTrait<T>
                                                             , DVecEvalExprTrait<T>
                                                             >::Type
                                                , typename If< IsTransposeVector<T>
                                                             , TSVecEvalExprTrait<T>
                                                             , SVecEvalExprTrait<T>
                                                             >::Type
                                                >::Type
                                   , Failure
                                   >::Type
                      >::Type  Tmp;

   typedef typename RemoveReference< typename RemoveCV<T>::Type >::Type  Type1;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, EvalExprTrait<Type1>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
