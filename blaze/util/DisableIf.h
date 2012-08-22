//=================================================================================================
/*!
//  \file blaze/util/DisableIf.h
//  \brief Header file for the DisableIf class template
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

#ifndef _BLAZE_UTIL_DISABLEIF_H_
#define _BLAZE_UTIL_DISABLEIF_H_


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
// The DisableIfTrue class template is an auxiliary tool for an intentional application of the
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
// every possible argument. However, with the DisableIfTrue class it is for example possible
// to prohibit built-in, numeric data types as argument types:

   \code
   template< typename Type >
   typename DisableIfTrue< IsNumeric<Type>::value >::Type process( Type t ) { ... }
   \endcode

// In case the given data type is a built-in, numeric data type, the access to the nested type
// definition \a Type of the DisableIfTrue class template will fail. However, due to the SFINAE
// principle, this will only result in a compilation error in case the compiler cannot find
// another valid function.\n
// Note that in this application of the DisableIfTrue template the default for the nested type
// definition \a Type is used, which corresponds to \a void. Via the second template argument
// it is possible to explicitly specify the type of \a Type:

   \code
   // Explicity specifying the default
   typename DisableIfTrue< IsNumeric<Type>::value, void >::Type

   // In case the given data type is not a boolean data type, the nested type definition
   // 'Type' is set to float
   typename DisableIfTrue< IsBoolean<Type>::value, float >::Type
   \endcode

// For more information on the DisableIfTrue/DisableIf functionality, see the Boost library
// documentation of the enable_if family at:
//
//           \a http://www.boost.org/doc/libs/1_40_0/libs/utility/enable_if.html.
*/
template< bool Condition     // Compile time condition
        , typename T=void >  // The type to be instantiated
struct DisableIfTrue
{
   //**********************************************************************************************
   typedef T  Type;  //!< The instantiated type.
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief DisableIfTrue specialization for failed constraints.
// \ingroup util
//
// This specialization of the DisableIfTrue template is selected if the first template parameter
// (the compile time condition) evaluates to \a true. This specialization does not contains a
// nested type definition \a Type and therefore always results in a compilation error in case
// \a Type is accessed. However, due to the SFINAE principle the compilation process is not
// necessarily stopped if another, valid instantiation is found by the compiler.
*/
template< typename T >  // The type to be instantiated
struct DisableIfTrue<true,T>
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
// The DisableIf class template is an auxiliary tool for an intentional application of the
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
// for every possible argument. However, with the DisableIf class it is for example possible
// to prohibit built-in, numeric data types as argument types:

   \code
   template< typename Type >
   typename DisableIf< IsNumeric<Type> >::Type process( Type t ) { ... }
   \endcode

// In case the given data type is a built-in, numeric data type, the access to the nested
// type definition \a Type of the DisableIf class template will fail. However, due to the
// SFINAE principle, this will only result in a compilation error in case the compiler cannot
// find another valid function.\n
// Note that in this application of the DisableIf template the default for the nested type
// definition \a Type is used, which corresponds to \a void. Via the second template argument
// it is possible to explicitly specify the type of \a Type:

   \code
   // Explicity specifying the default
   typename DisableIf< IsNumeric<Type>, void >::Type

   // In case the given data type is not a boolean data type, the nested type definition
   // 'Type' is set to float
   typename DisableIf< IsBoolean<Type>, float >::Type
   \endcode

// Note that in contrast to the DisableIfTrue template, the DisableIf template expects a
// type as first template argument that has a nested type definition \a value. Therefore
// the DisableIf template is the more convenient choice for all kinds of type traits.
//
// For more information on the DisableIfTrue/DisableIf functionality, see the Boost library
// documentation of the enable_if family at:
//
//           \a http://www.boost.org/doc/libs/1_40_0/libs/utility/enable_if.html.
*/
template< typename Condition  // Compile time condition
        , typename T=void >   // The type to be instantiated
struct DisableIf : public DisableIfTrue<Condition::value,T>
{};
//*************************************************************************************************

} // namespace blaze

#endif
