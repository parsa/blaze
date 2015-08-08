//=================================================================================================
/*!
//  \file blaze/util/StaticAssert.h
//  \brief Compile time assertion
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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

// \note: blaze::BLAZE_STATIC_ASSERT can only can expressions of integral type. Floating point
// expressions cannot be checked at compile time!
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
