//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsVectorizable.h
//  \brief Header file for the IsVectorizable type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISVECTORIZABLE_H_
#define _BLAZE_UTIL_TYPETRAITS_ISVECTORIZABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Vectorization.h>
#include <blaze/util/Complex.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsVectorizable type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsVectorizableHelper
{
   //**********************************************************************************************
   enum { value = ( BLAZE_SSE_MODE  && ( IsFloat<T>::value  || IsSame<complex<float>,T>::value  ) ) ||
                  ( BLAZE_SSE2_MODE && ( IsDouble<T>::value || IsSame<complex<double>,T>::value ) ) ||
                  ( BLAZE_SSE2_MODE && ( IsIntegral<T>::value && sizeof(T) >= 2UL ) ) ||
                  ( BLAZE_MIC_MODE  && ( IsIntegral<T>::value && sizeof(T) >= 4UL ) ) ||
                  ( BLAZE_MIC_MODE  && ( IsFloat<T>::value  || IsSame<complex<float>,T>::value  ) ) ||
                  ( BLAZE_MIC_MODE  && ( IsDouble<T>::value || IsSame<complex<double>,T>::value ) ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for vectorizable types.
// \ingroup type_traits
//
// Depending on the available instruction set (SSE, SSE2, SSE3, SSE4, AVX, AVX2, MIC, ...),
// this type trait tests whether or not the given template parameter is a vectorizable type,
// i.e. a type for which intrinsic vector operations and optimizations can be used. Currently,
// only signed/unsigned short, signed/unsigned int, signed/unsigned long, float, double,
// complex<float>, and complex<double> are considered to be vectorizable types. In case the
// type is vectorizable, the \a value member enumeration is set to 1, the nested type definition
// \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to
// 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsVectorizable< int >::value         // Evaluates to 1
   blaze::IsVectorizable< const float >::Type  // Results in TrueType
   blaze::IsVectorizable< volatile double >    // Is derived from TrueType
   blaze::IsVectorizable< bool >::value        // Evaluates to 0
   blaze::IsVectorizable< const char >::Type   // Results in FalseType
   blaze::IsVectorizable< long double >        // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsVectorizable : public IsVectorizableHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsVectorizableHelper<T>::value };
   typedef typename IsVectorizableHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
