//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsVector.h
//  \brief Header file for the IsVector type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISVECTOR_H_
#define _BLAZE_MATH_TYPETRAITS_ISVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
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
/*!\brief Auxiliary helper struct for the IsVector type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsVectorHelper
{
   //**********************************************************************************************
   enum { value = IsDenseVector<T>::value || IsSparseVector<T>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for vector types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a N-dimensional dense
// or sparse vector type. In case the type is a vector type, the \a value member enumeration
// is set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::IsVector< StaticVector<float,3U,false> >::value      // Evaluates to 1
   blaze::IsVector< const DynamicVector<double,true> >::Type   // Results in TrueType
   blaze::IsVector< volatile CompressedVector<int,true> >      // Is derived from TrueType
   blaze::IsVector< StaticMatrix<double,3U,3U,false> >::value  // Evaluates to 0
   blaze::IsVector< const DynamicMatrix<double,true> >::Type   // Results in FalseType
   blaze::IsVector< volatile CompressedMatrix<int,true> >      // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsVector : public IsVectorHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsVectorHelper<T>::value };
   typedef typename IsVectorHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
