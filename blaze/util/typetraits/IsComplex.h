//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsComplex.h
//  \brief Header file for the IsComplex type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISCOMPLEX_H_
#define _BLAZE_UTIL_TYPETRAITS_ISCOMPLEX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for complex types.
// \ingroup type_traits
//
// This type trait tests whether or not the given template parameter is a complex data type.
// In case the type is a complex data type (ignoring the cv-qualifiers), the \a value member
// enumeration is set to 1, the nested type definition \a Type is \a TrueType, and the class
// derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the
// class derives from \a FalseType.

   \code
   blaze::IsComplex< complex<double> >::value      // Evaluates to 1
   blaze::IsComplex< const complex<float> >::Type  // Results in TrueType
   blaze::IsComplex< volatile complex<int> >       // Is derived from TrueType
   blaze::IsComplex< float >::value                // Evaluates to 0
   blaze::IsComplex< const double >::Type          // Results in FalseType
   blaze::IsComplex< const volatile int >          // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsComplex : public FalseType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = 0 };
   typedef FalseType  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsComplex type trait for the plain 'complex' type.
template< typename T >
struct IsComplex< complex<T> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsComplex type trait for 'const complex'.
template< typename T >
struct IsComplex< const complex<T> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsComplex type trait for 'volatile complex'.
template< typename T >
struct IsComplex< volatile complex<T> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsComplex type trait for 'const volatile complex'
template< typename T >
struct IsComplex< const volatile complex<T> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
