//=================================================================================================
/*!
//  \file blaze/util/typetraits/HaveSameSize.h
//  \brief Header file for the HaveSameSize type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_HAVESAMESIZE_H_
#define _BLAZE_UTIL_TYPETRAITS_HAVESAMESIZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

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
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This class offers the possibility to test the size of two types at compile time. If an object
// of type \a T1 has the same size as an object of type \a T2, the \a value member enumeration is
// set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value  is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::HaveSameSize<int,unsigned int>::value  // Evaluates to 1
   blaze::HaveSameSize<int,unsigned int>::Type   // Results in TrueType
   blaze::HaveSameSize<int,unsigned int>         // Is derived from TrueType
   blaze::HaveSameSize<char,wchar_t>::value      // Evalutes to 0
   blaze::HaveSameSize<char,wchar_t>::Type       // Results in FalseType
   blaze::HaveSameSize<char,wchar_t>             // Is derived from FalseType
   \endcode

// One example for the application of this type trait is a compile time check if the compiler
// supports the 'Empty Derived class Optimization (EDO)':

   \code
   // Definition of the base class A
   struct A {
      int i_;
   };

   // Definition of the derived class B
   struct B : public A {};

   // Testing whether or not an object of type B has the same size as the
   //   base class A and whether the compiler supports EDO
   blaze::HaveSameSize( A, B );
   \endcode
*/
template< typename T1, typename T2 >
class HaveSameSize : public SelectType< sizeof(T1) == sizeof(T2), TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = sizeof( T1 ) == sizeof( T2 ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// // This class is a partial specialization of the HaveSameSize template for the type \a void
// as first template argument. The \a value member enumeration is automatically set to 0, the
// nested type definition \a Type is \a FalseType, and the class derives from \a FalseType for
// any given type \a T since the \a void type has no size.
*/
template< typename T >
class HaveSameSize<void,T> : public FalseType
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class is a partial specialization of the HaveSameSize template for the type \a void
// as second template argument. The \a value member enumeration is automatically set to 0, the
// nested type definition \a Type is \a FalseType, and the class derives from \a FalseType for
// any given type \a T since the \a void type has no size.
*/
template< typename T >
class HaveSameSize<T,void> : public FalseType
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Full specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class is a full specialization of the HaveSameSize template for the type \a void
// as first and second template argument. The \a value member enumeration is automatically
// set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType since both arguments are \a void.
*/
template<>
class HaveSameSize<void,void> : public TrueType
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
