//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsMultExpr.h
//  \brief Header file for the IsMultExpr type trait class
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISMULTEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_base_of.hpp>
#include <blaze/math/expressions/MultExpr.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsMultExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsMultExprHelper
{
   //**********************************************************************************************
   enum { value = boost::is_base_of<MultExpr,T>::value && !boost::is_base_of<T,MultExpr>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check whether the given type is a multiplication expression template.
// \ingroup math_type_traits
//
// This type trait class tests whether or not the given type \a Type is a multiplication
// expression template (i.e. an expression representing an element-wise vector multiplication,
// a matrix/vector multiplication, a vector/matrix multiplication, or a matrix multiplication).
// In order to qualify as a valid multiplication expression template, the given type has to
// derive (publicly or privately) from the MultExpr base class. In case the given type is a
// valid multiplication expression template, the \a value member enumeration is set to 1, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.
*/
template< typename T >
struct IsMultExpr : public IsMultExprHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsMultExprHelper<T>::value };
   typedef typename IsMultExprHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
