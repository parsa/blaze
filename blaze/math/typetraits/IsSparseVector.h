//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsSparseVector.h
//  \brief Header file for the IsSparseVector type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSPARSEVECTOR_H_
#define _BLAZE_MATH_TYPETRAITS_ISSPARSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_base_of.hpp>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsSparseVector type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSparseVectorHelper
{
 private:
   //**********************************************************************************************
   typedef typename RemoveCV<T>::Type  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { value = boost::is_base_of< SparseVector<T2,false>, T2 >::value ||
                  boost::is_base_of< SparseVector<T2,true >, T2 >::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for sparse vector types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a sparse, N-dimensional
// vector type. In case the type is a sparse vector type, the \a value member enumeration is
// set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::IsSparseVector< CompressedVector<double,false> >::value     // Evaluates to 1
   blaze::IsSparseVector< const CompressedVector<float,true> >::Type  // Results in TrueType
   blaze::IsSparseVector< volatile CompressedVector<int,true> >       // Is derived from TrueType
   blaze::IsSparseVector< DynamicVector<double,false> >::value        // Evaluates to 0
   blaze::IsSparseVector< const DynamicMatrix<double,true> >::Type    // Results in FalseType
   blaze::IsSparseVector< CompressedMatrix<double,true> >             // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsSparseVector : public IsSparseVectorHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsSparseVectorHelper<T>::value };
   typedef typename IsSparseVectorHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
