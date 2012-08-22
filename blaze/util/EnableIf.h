//=================================================================================================
/*!
//  \file blaze/util/EnableIf.h
//  \brief Header file for the EnableIf class template
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

#ifndef _BLAZE_UTIL_ENABLEIF_H_
#define _BLAZE_UTIL_ENABLEIF_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Substitution Failure Is Not An Error (SFINAE) class.
// \ingroup util
//
// The EnableIfTrue class template is an auxiliary tool for an intentional application of the
// Substitution Failure Is Not An Error (SFINAE) principle. It allows a function template or a
// class template specialization to include or exclude itself from a set of matching functions
// or specializations based on properties of its template arguments. For instance, it can be
// used to restrict the selection of a function template to specific data types. The following
// example illustrates this in more detail.

   \code
   template< typename Type >
   void process( Type t ) { ... }
   \endcode

// Due to the general formulation of this function, it will always be a possible candidate for
// every possible argument. However, with the EnableIfTrue class it is for example possible to
// restrict the valid argument types to built-in, numeric data types.

   \code
   template< typename Type >
   typename EnableIfTrue< IsNumeric<Type>::value >::Type process( Type t ) { ... }
   \endcode

// In case the given data type is not a built-in, numeric data type, the access to the nested
// type defintion \a Type of the EnableIfTrue template will fail. However, due to the SFINAE
// principle, this will only result in a compilation error in case the compiler cannot find
// another valid function.\n
// Note that in this application of the EnableIfTrue template the default for the nested type
// definition \a Type is used, which corresponds to \a void. Via the second template argument
// it is possible to explicitly specify the type of \a Type:

   \code
   // Explicity specifying the default
   typename EnableIfTrue< IsNumeric<Type>::value, void >::Type

   // In case the given data type is a boolean data type, the nested type definition
   // 'Type' is set to float
   typename EnableIfTrue< IsBoolean<Type>::value, float >::Type
   \endcode

// For more information on the EnableIfTrue/EnableIf functionality, see the Boost library
// documentation of the enable_if family at:
//
//           \a http://www.boost.org/doc/libs/1_40_0/libs/utility/enable_if.html.
*/
template< bool Condition     // Compile time condition
        , typename T=void >  // The type to be instantiated
struct EnableIfTrue
{
   //**********************************************************************************************
   typedef T  Type;  //!< The instantiated type.
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief EnableIfTrue specialization for failed constraints.
// \ingroup util
//
// This specialization of the EnableIfTrue template is selected if the first template parameter
// (the compile time condition) evaluates to \a false. This specialization does not contains a
// nested type definition \a Type and therefore always results in a compilation error in case
// \a Type is accessed. However, due to the SFINAE principle the compilation process is not
// necessarily stopped if another, valid instantiation is found by the compiler.
*/
template< typename T >  // The type to be instantiated
struct EnableIfTrue<false,T>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Substitution Failure Is Not An Error (SFINAE) class.
// \ingroup util
//
// The EnableIf class template is an auxiliary tool for an intentional application of the
// Substitution Failure Is Not An Error (SFINAE) principle. It allows a function template
// or a class template specialization to include or exclude itself from a set of matching
// functions or specializations based on properties of its template arguments. For instance,
// it can be used to restrict the selection of a function template to specific data types.
// The following example illustrates this in more detail.

   \code
   template< typename Type >
   void process( Type t ) { ... }
   \endcode

// Due to the general formulation of this function, it will always be a possible candidate
// for every possible argument. However, with the EnableIf class it is for example possible
// to restrict the valid argument types to built-in, numeric data types.

   \code
   template< typename Type >
   typename EnableIf< IsNumeric<Type> >::Type process( Type t ) { ... }
   \endcode

// In case the given data type is not a built-in, numeric data type, the access to the nested
// type defintion \a Type of the EnableIf template will fail. However, due to the SFINAE
// principle, this will only result in a compilation error in case the compiler cannot find
// another valid function.\n
// Note that in this application of the EnableIf template the default for the nested type
// definition \a Type is used, which corresponds to \a void. Via the second template argument
// it is possible to explicitly specify the type of \a Type:

   \code
   // Explicity specifying the default
   typename EnableIf< IsNumeric<Type>, void >::Type

   // In case the given data type is a boolean data type, the nested type definition
   // 'Type' is set to float
   typename EnableIf< IsBoolean<Type>, float >::Type
   \endcode

// Note that in contrast to the EnableIfTrue template, the EnableIf template expects a type as
// first template argument that has a nested type definition \a value. Therefore the EnableIf
// template is the more convenient choice for all kinds of type traits.
//
// For more information on the EnableIfTrue/EnableIf functionality, see the Boost library
// documentation of the enable_if family at:
//
//           \a http://www.boost.org/doc/libs/1_40_0/libs/utility/enable_if.html.
*/
template< typename Condition     // Compile time condition
        , typename T=void >  // The type to be instantiated
struct EnableIf : public EnableIfTrue<Condition::value,T>
{};
//*************************************************************************************************

} // namespace blaze

#endif
