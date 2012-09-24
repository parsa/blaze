//=================================================================================================
/*!
//  \file blaze/math/traits/TDMatScalarDivTrait.h
//  \brief Header file for the TDMatScalarDivTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_TDMATSCALARDIVTRAIT_H_
#define _BLAZE_MATH_TRAITS_TDMATSCALARDIVTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the TDMatScalarDivTrait trait.
// \ingroup math_traits
*/
template< typename MT
        , typename ST
        , bool Condition >
struct TDMatScalarDivTraitHelper
{
 private:
   //**********************************************************************************************
   typedef typename DivTrait<typename BaseElementType<MT>::Type,ST>::Type  ElementType;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsFloatingPoint<ElementType>::value
                              , DMatScalarMultExpr<MT,ElementType,true>
                              , DMatScalarDivExpr<MT,ElementType,true> >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TDMatScalarDivTraitHelper class template.
// \ingroup math_traits
*/
template< typename MT
        , typename ST >
struct TDMatScalarDivTraitHelper<MT,ST,false>
{
 public:
   //**********************************************************************************************
   typedef INVALID_TYPE  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluation of the expression type of a transpose dense matrix/scalar division.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a transpose
// dense matrix/scalar division. Given the column-major dense matrix type \a MT and the scalar
// type \a ST, the nested type \a Type corresponds to the resulting expression type. In case
// either \a MT is not a column-major dense matrix type or \a ST is not a scalar type, the
// resulting \a Type is set to \a INVALID_TYPE.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename ST >  // Type of the right-hand side scalar
struct TDMatScalarDivTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<MT>::value || IsVolatile<MT>::value || IsReference<MT>::value ||
                      IsConst<ST>::value || IsVolatile<ST>::value || IsReference<ST>::value };

   enum { condition = IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                      IsNumeric<ST>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef TDMatScalarDivTraitHelper<MT,ST,condition>  Tmp;

   typedef typename RemoveReference< typename RemoveCV<MT>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<ST>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, TDMatScalarDivTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
