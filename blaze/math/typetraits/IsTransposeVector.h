//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsTransposeVector.h
//  \brief Header file for the IsTransposeVector type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISTRANSPOSEVECTOR_H_
#define _BLAZE_MATH_TYPETRAITS_ISTRANSPOSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename, bool > struct DenseVector;
template< typename, bool > struct SparseVector;




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsTransposeVector type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsTransposeVectorHelper
{
 private:
   //**********************************************************************************************
   typedef typename boost::remove_cv<T>::type  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { value = boost::is_base_of< DenseVector <T2,true>, T2 >::value ||
                  boost::is_base_of< SparseVector<T2,true>, T2 >::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for transpose vector types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template argument is a transpose dense or
// sparse vector type (i.e., a vector whose transposition flag is set to 1). In case the type
// is a transpose vector type, the \a value member enumeration is set to 1, the nested type
// definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsTransposeVector< StaticVector<float,3U,true> >::value       // Evaluates to 1
   blaze::IsTransposeVector< const DynamicVector<double,true> >::Type   // Results in TrueType
   blaze::IsTransposeVector< volatile CompressedVector<int,true> >      // Is derived from TrueType
   blaze::IsTransposeVector< StaticVector<float,3U,false> >::value      // Evaluates to 0
   blaze::IsTransposeVector< const DynamicVector<double,false> >::Type  // Results in FalseType
   blaze::IsTransposeVector< volatile CompressedVector<int,false> >     // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsTransposeVector : public IsTransposeVectorHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsTransposeVectorHelper<T>::value };
   typedef typename IsTransposeVectorHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
