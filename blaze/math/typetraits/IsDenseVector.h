//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsDenseVector.h
//  \brief Header file for the IsDenseVector type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISDENSEVECTOR_H_
#define _BLAZE_MATH_TYPETRAITS_ISDENSEVECTOR_H_


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




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsDenseVector type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsDenseVectorHelper
{
 private:
   //**********************************************************************************************
   typedef typename boost::remove_cv<T>::type  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { value = boost::is_base_of< DenseVector<T2,false>, T2 >::value ||
                  boost::is_base_of< DenseVector<T2,true >, T2 >::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for dense vector types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a dense, N-dimensional
// vector type. In case the type is a dense vector type, the \a value member enumeration is
// set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::IsDenseVector< DynamicVector<double,false> >::value       // Evaluates to 1
   blaze::IsDenseVector< const StaticVector<float,3U,true> >::Type  // Results in TrueType
   blaze::IsDenseVector< volatile StaticVector<int,6U,true> >       // Is derived from TrueType
   blaze::IsDenseVector< CompressedVector<double,false> >::value    // Evaluates to 0
   blaze::IsDenseVector< CompressedMatrix<double,true> >::Type      // Results in FalseType
   blaze::IsDenseVector< DynamicMatrix<double,true> >               // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsDenseVector : public IsDenseVectorHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsDenseVectorHelper<T>::value };
   typedef typename IsDenseVectorHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
