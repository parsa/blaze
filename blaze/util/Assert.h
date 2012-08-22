//=================================================================================================
/*!
//  \file blaze/util/Assert.h
//  \brief Header file for run time assertion macros
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

#ifndef _BLAZE_UTIL_ASSERT_H_
#define _BLAZE_UTIL_ASSERT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cassert>
#include <blaze/system/Assertion.h>


namespace blaze {

//=================================================================================================
//
//  RUN TIME ASSERTION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup assert Assertions
// \ingroup util
*/
/*!\defgroup runtime_assert Run time assertions
// \ingroup assert
*/
/*!\brief Assertion helper function.
// \ingroup runtime_assert
//
// The ASSERT_MESSAGE function is a small helper function to assist in printing an informative
// message in case an assert fires. This function builds on the ideas of Matthew Wilson, who
// directly combines a C-string error message with the run time expression (Imperfect C++,
// ISBN: 0321228774):

   \code
   assert( ... &&  "Error message" );
   assert( ... || !"Error message" );
   \endcode

// However, both approaches fail to compile without warning on certain compilers. Therefore
// this inline function is used instead of the direct approaches, which circumvents all compiler
// warnings:

   \code
   assert( ... || ASSERT_MESSAGE( "Error message" ) );
   \endcode
*/
inline bool ASSERT_MESSAGE( const char* /*msg*/ )
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Run time assertion macro for internal checks.
// \ingroup runtime_assert
//
// In case of an invalid run time expression, the program execution is terminated.\n
// The BLAZE_INTERNAL_ASSERT macro can be disabled by setting the \a BLAZE_USER_ASSERTION
// flag to zero or by defining \a NDEBUG during the compilation.
*/
#if BLAZE_INTERNAL_ASSERTION
#  define BLAZE_INTERNAL_ASSERT(expr,msg) assert( ( expr ) || blaze::ASSERT_MESSAGE( msg ) )
#else
#  define BLAZE_INTERNAL_ASSERT(expr,msg)
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Run time assertion macro for user checks.
// \ingroup runtime_assert
//
// In case of an invalid run time expression, the program execution is terminated.\n
// The BLAZE_USER_ASSERT macro can be disabled by setting the \a BLAZE_USER_ASSERT flag
// to zero or by defining \a NDEBUG during the compilation.
*/
#if BLAZE_USER_ASSERTION
#  define BLAZE_USER_ASSERT(expr,msg) assert( ( expr ) || blaze::ASSERT_MESSAGE( msg ) )
#else
#  define BLAZE_USER_ASSERT(expr,msg)
#endif
//*************************************************************************************************

} // namespace blaze

#endif
