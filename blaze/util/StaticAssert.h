//=================================================================================================
/*!
//  \file blaze/util/StaticAssert.h
//  \brief Compile time assertion
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

#ifndef _BLAZE_UTIL_STATICASSERT_H_
#define _BLAZE_UTIL_STATICASSERT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Suffix.h>


namespace blaze {

//=================================================================================================
//
//  COMPILE TIME ASSERTION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup static_assert Compile time assertion
// \ingroup assert
//
// Static assertion offers the possibility to stop the compilation process if a specific
// compile time condition is not met. The blaze::BLAZE_STATIC_ASSERT macro can be used to check
// an integral constant expression at compile time. If the expression evaluates to \a false,
// a compilation error is generated that stops the compilation process. If the expression
// (hopefully) evaluates to \a true, the compilation process is not aborted and the static
// check leaves neither code nor data and is therefore not affecting the performance.\n
// The blaze::BLAZE_STATIC_ASSERT macro can be used wherever a standard typedef statement can
// be declared, i.e. in namespace scope, in class scope and in function scope. The following
// examples illustrate the use of the static assertion macro: the type of the rotation matrix
// is checked at compile time and restricted to be of floating point type.

   \code
   #include <limits>

   template< typename T >
   class RotationMatrix {
      ...
      BLAZE_STATIC_ASSERT( !std::numeric_limits<T>::is_integer );
      ...
   };
   \endcode

// The static assertion is implemented in such a way that the created error messages for a
// failed compile time check contains either the term STATIC_ASSERT or STATIC_ASSERT_FAILED.
// The error message doesn't explicitly explain the source of the error, but is at least
// useful to catch the eye. The following examples show possible error messages for the
// GNU g++ and the Intel compiler:

   \code
   incomplete type ‘blaze::STATIC_ASSERTION_FAILED<false>’ used in nested name specifier
   \endcode

   \code
   error: incomplete type is not allowed
      BLAZE_STATIC_ASSERT( !std::numeric_limits<T>::is_integer );
   \endcode

// \b Note: blaze::BLAZE_STATIC_ASSERT can only can expressions of integral type. Floating
// point expressions cannot be checked at compile time!
//
// \b Acknowledgements: blaze::BLAZE_STATIC_ASSERT builds on ideas developed by John Maddock
// within the Boost C++ framework (www.boost.org). However, it uses a slightly changed compile
// time check that completely relies on a nested enum variable and doesn't use an old-style C
// cast.
*/
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Wrapper class for the static assertion check.
// \ingroup static_assert
//
// This class is used as a wrapper for the instantiation of the STATIC_ASSERTION_FAILED
// class template. It serves the purpose to force the instantiation of either the defined
// specialization or the undefined basic template during the compilation. In case the
// compile time condition is met, the type blaze::STATIC_ASSERTION_TEST<1> is defined.
*/
template< int > struct STATIC_ASSERTION_TEST {};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time assertion
// \ingroup static_assert
//
// Helper template class for the compile time assertion. Based on the compile time constant
// expression used for the template instantiation, either the undefined basic template or the
// specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct STATIC_ASSERTION_FAILED;
template<> struct STATIC_ASSERTION_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time assertion macro.
// \ingroup static_assert
//
// In case of an invalid compile time expression, a compilation error is created.
*/
#define BLAZE_STATIC_ASSERT(expr) \
   typedef ::blaze::STATIC_ASSERTION_TEST< ::blaze::STATIC_ASSERTION_FAILED< (expr) != 0 >::value > \
      BLAZE_JOIN( BLAZE_STATIC_ASSERTION_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
