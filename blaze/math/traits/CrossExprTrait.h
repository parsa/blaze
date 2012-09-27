//=================================================================================================
/*!
//  \file blaze/math/traits/CrossExprTrait.h
//  \brief Header file for the CrossExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_CROSSEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_CROSSEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DVecDVecCrossExprTrait.h>
#include <blaze/math/traits/DVecSVecCrossExprTrait.h>
#include <blaze/math/traits/SVecDVecCrossExprTrait.h>
#include <blaze/math/traits/SVecSVecCrossExprTrait.h>
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
/*!\brief Evaluation of the return type of a cross product expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a cross product expression.
// Given the two types \a T1 and \a T2, which must be dense or sparse column vectors, the nested
// type \a Type corresponds to the resulting return type. In case the types of \a T1 or \a T2
// don't fit or if the two types cannot be used in a cross product, the resulting data type
// \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side cross product operand
        , typename T2 >  // Type of the right-hand side cross product operand
struct CrossExprTrait
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
   typedef typename If< IsVector<T1>
                      , typename If< IsVector<T2>
                                   , typename IfNot< IsTransposeVector<T1>
                                                   , typename IfNot< IsTransposeVector<T2>
                                                                   , typename If< IsDenseVector<T1>
                                                                                , typename If< IsDenseVector<T2>
                                                                                             , DVecDVecCrossExprTrait<T1,T2>
                                                                                             , DVecSVecCrossExprTrait<T1,T2>
                                                                                             >::Type
                                                                                , typename If< IsDenseVector<T2>
                                                                                             , SVecDVecCrossExprTrait<T1,T2>
                                                                                             , SVecSVecCrossExprTrait<T1,T2>
                                                                                             >::Type
                                                                                >::Type
                                                                   , Failure
                                                                   >::Type
                                                   , Failure
                                                   >::Type
                                   , Failure
                                   >::Type
                      , Failure
                      >::Type  Tmp;

   typedef typename RemoveReference< typename RemoveCV<T1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<T2>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, CrossExprTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
