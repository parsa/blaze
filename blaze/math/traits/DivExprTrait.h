//=================================================================================================
/*!
//  \file blaze/math/traits/DivExprTrait.h
//  \brief Header file for the DivExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_DIVEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_DIVEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/DMatScalarDivExprTrait.h>
#include <blaze/math/traits/DVecScalarDivExprTrait.h>
#include <blaze/math/traits/SMatScalarDivExprTrait.h>
#include <blaze/math/traits/SVecScalarDivExprTrait.h>
#include <blaze/math/traits/TDMatScalarDivExprTrait.h>
#include <blaze/math/traits/TDVecScalarDivExprTrait.h>
#include <blaze/math/traits/TSMatScalarDivExprTrait.h>
#include <blaze/math/traits/TSVecScalarDivExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
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
/*!\brief Evaluation of the resulting expression type of a division.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a division expression
// between scalars, vectors, and matrices. Given the two types \a T1 and \a T2, where \a T1 must
// either be a scalar, vector, or matrix type and \a T2 which must be a scalar type, the nested
// type \a Type corresponds to the resulting return type. In case \a T1 or \a T2 don't fit or if
// the two types cannot be divided, the resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side division operand
        , typename T2 >  // Type of the right-hand side division operand
struct DivExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<T1>::value || IsVolatile<T1>::value || IsReference<T1>::value ||
                      IsConst<T2>::value || IsVolatile<T2>::value || IsReference<T2>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsMatrix<T1>
                      , typename If< IsDenseMatrix<T1>
                                   , typename If< IsRowMajorMatrix<T1>
                                                , typename If< IsNumeric<T2>
                                                             , DMatScalarDivExprTrait<T1,T2>
                                                             , Failure
                                                             >::Type
                                                , typename If< IsNumeric<T2>
                                                             , TDMatScalarDivExprTrait<T1,T2>
                                                             , Failure
                                                             >::Type
                                                >::Type
                                   , typename If< IsRowMajorMatrix<T1>
                                                , typename If< IsNumeric<T2>
                                                             , SMatScalarDivExprTrait<T1,T2>
                                                             , Failure
                                                             >::Type
                                                , typename If< IsNumeric<T2>
                                                             , TSMatScalarDivExprTrait<T1,T2>
                                                             , Failure
                                                             >::Type
                                                >::Type
                                   >::Type
                      , typename If< IsVector<T1>
                                   , typename If< IsDenseVector<T1>
                                                , typename If< IsTransposeVector<T1>
                                                             , typename If< IsNumeric<T2>
                                                                          , TDVecScalarDivExprTrait<T1,T2>
                                                                          , Failure
                                                                          >::Type
                                                             , typename If< IsNumeric<T2>
                                                                          , DVecScalarDivExprTrait<T1,T2>
                                                                          , Failure
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsTransposeVector<T1>
                                                             , typename If< IsNumeric<T2>
                                                                          , TSVecScalarDivExprTrait<T1,T2>
                                                                          , Failure
                                                                          >::Type
                                                             , typename If< IsNumeric<T2>
                                                                          , SVecScalarDivExprTrait<T1,T2>
                                                                          , Failure
                                                                          >::Type
                                                             >::Type
                                                >::Type
                                   , typename If< IsNumeric<T1>
                                                , typename If< IsNumeric<T2>
                                                             , DivTrait<T1,T2>
                                                             , Failure
                                                             >::Type
                                                , Failure
                                                >::Type
                                   >::Type
                      >::Type  Tmp;

   typedef typename RemoveReference< typename RemoveCV<T1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<T2>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, DivExprTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
