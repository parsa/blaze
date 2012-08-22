//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsBlasCompatible.h
//  \brief Header file for the IsBlasCompatible type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISBLASCOMPATIBLE_H_
#define _BLAZE_MATH_TYPETRAITS_ISBLASCOMPATIBLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsBlasCompatible type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsBlasCompatibleHelper
{
   //**********************************************************************************************
   enum { value = ( IsFloat<T>::value  || IsSame<complex<float>,T>::value ||
                    IsDouble<T>::value || IsSame<complex<double>,T>::value ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for data types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a data type compatible
// to the BLAS standard. The BLAS standard currently only supports float, double, complex<float>
// and complex<double> values. If the type is BLAS compatible, the \a value member enumeration
// is set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::IsBlasCompatible< float >::value         // Evaluates to 1
   blaze::IsBlasCompatible< double >::Type         // Results in TrueType
   blaze::IsBlasCompatible< complex<float> >       // Is derived from TrueType
   blaze::IsBlasCompatible< int >::value           // Evaluates to 0
   blaze::IsBlasCompatible< unsigned long >::Type  // Results in FalseType
   blaze::IsBlasCompatible< long double >          // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsBlasCompatible : public IsBlasCompatibleHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsBlasCompatibleHelper<T>::value };
   typedef typename IsBlasCompatibleHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
