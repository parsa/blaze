//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsNumeric.h
//  \brief Header file for the IsNumeric type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISNUMERIC_H_
#define _BLAZE_UTIL_TYPETRAITS_ISNUMERIC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/typetraits/IsBoolean.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsVoid.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsNumeric type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsNumericHelper
{
   //**********************************************************************************************
   enum { value = ( IsBuiltin<T>::value && !IsBoolean<T>::value && !IsVoid<T>::value ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for numeric types.
// \ingroup type_traits
//
// This type trait tests whether or not the given template parameter is a numeric data type.
// Blaze considers all integral (except \a bool), floating point, and complex data types as
// numeric data types. In case the type is a numeric type, the \a value member enumeration is
// set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::IsNumeric<int>::value                // Evaluates to 1 (int is a numeric data type)
   blaze::IsNumeric<const double>::Type        // Results in TrueType (float is a numeric data type)
   blaze::IsNumeric<volatile complex<float> >  // Is derived from TrueType (complex<float> is a numeric data type)
   blaze::IsNumeric<void>::value               // Evaluates to 0 (void is not a numeric data type)
   blaze::IsNumeric<bool>::Type                // Results in FalseType (bool is not a numeric data type)
   blaze::IsNumeric<const bool>                // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsNumeric : public IsNumericHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsNumericHelper<T>::value };
   typedef typename IsNumericHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsNumeric type trait for the plain 'complex' type.
template< typename T >
struct IsNumeric< complex<T> > : public IsNumeric<T>::Type
{
 public:
   //**********************************************************************************************
   enum { value = IsNumeric<T>::value };
   typedef typename IsNumeric<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsNumeric type trait for 'const complex'.
template< typename T >
struct IsNumeric< const complex<T> > : public IsNumeric<T>::Type
{
 public:
   //**********************************************************************************************
   enum { value = IsNumeric<T>::value };
   typedef typename IsNumeric<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsNumeric type trait for 'volatile complex'.
template< typename T >
struct IsNumeric< volatile complex<T> > : public IsNumeric<T>::Type
{
 public:
   //**********************************************************************************************
   enum { value = IsNumeric<T>::value };
   typedef typename IsNumeric<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsNumeric type trait for 'const volatile complex'.
template< typename T >
struct IsNumeric< const volatile complex<T> > : public IsNumeric<T>::Type
{
 public:
   //**********************************************************************************************
   enum { value = IsNumeric<T>::value };
   typedef typename IsNumeric<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
