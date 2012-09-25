//=================================================================================================
/*!
//  \file blaze/math/traits/TSVecDVecMultExprTrait.h
//  \brief Header file for the TSVecDVecMultExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_TSVECDVECMULTEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_TSVECDVECMULTEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the TSVecDVecMultExprTrait.
// \ingroup math_traits
*/
template< typename VT1
        , typename VT2
        , bool Valid >
struct TSVecDVecMultExprTraitHelper
{
   //**********************************************************************************************
   typedef INVALID_TYPE  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TSVecDVecMultExprTrait auxiliary helper struct.
// \ingroup math_traits
*/
template< typename VT1
        , typename VT2 >
struct TSVecDVecMultExprTraitHelper<VT1,VT2,true>
{
   //**********************************************************************************************
   typedef typename MultTrait<typename VT1::ElementType,typename VT2::ElementType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluation of the expression type of a transpose sparse vector/dense vector multiplication.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a transpose
// sparse vector/dense vector multiplication (inner product). Given the transpose sparse vector
// type \a VT1 and the non-transpose dense vector type \a VT2, the nested type \a Type corresponds
// to the resulting expression type. In case either \a VT1 is not a transpose sparse vector
// type or \a VT2 is not a non-transpose dense vector type, the resulting \a Type is set to
// \a INVALID_TYPE.
*/
template< typename VT1    // Type of the left-hand side transpose sparse vector
        , typename VT2 >  // Type of the right-hand side non-transpose dense vector
struct TSVecDVecMultExprTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<VT1>::value || IsVolatile<VT1>::value || IsReference<VT1>::value ||
                      IsConst<VT2>::value || IsVolatile<VT2>::value || IsReference<VT2>::value };

   enum { valid = IsSparseVector<VT1>::value && IsTransposeVector<VT1>::value &&
                  IsDenseVector<VT2>::value  && !IsTransposeVector<VT2>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef TSVecDVecMultExprTraitHelper<VT1,VT2,valid>  Tmp;

   typedef typename RemoveReference< typename RemoveCV<VT1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<VT2>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, TSVecDVecMultExprTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
