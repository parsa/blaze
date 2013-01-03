//=================================================================================================
/*!
//  \file blaze/math/typetraits/RequiresEvaluation.h
//  \brief Header file for the RequiresEvaluation type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_REQUIRESEVALUATION_H_
#define _BLAZE_MATH_TYPETRAITS_REQUIRESEVALUATION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the RequiresEvaluation type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct RequiresEvaluationHelper
{
   //**********************************************************************************************
   enum { value = !IsReference<typename T::CompositeType>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check to query the requirement to evaluate an expression.
// \ingroup math_type_traits
//
// Via this type trait it is possible to determine whether a given vector or matrix expression
// type requires an intermediate evaluation in the context of a compound expression. In case
// the given type requires an evaluation, the \a value member enumeration is set to 1, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.
//
// \b Note that this type trait can only be applied to Blaze vector or matrix expressions
// or any other type providing the nested type \a CompositeType. In case this nested type
// is not available, applying the type trait results in a compile time error!
*/
template< typename T >
struct RequiresEvaluation : public RequiresEvaluationHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = RequiresEvaluationHelper<T>::value };
   typedef typename RequiresEvaluationHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
