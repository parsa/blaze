//=================================================================================================
/*!
//  \file blaze/util/constraints/SameType.h
//  \brief Data type constraint
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

#ifndef _BLAZE_UTIL_CONSTRAINTS_SAMETYPE_H_
#define _BLAZE_UTIL_CONSTRAINTS_SAMETYPE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/ConstraintTest.h>
#include <blaze/util/Suffix.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  MUST_BE_SAME_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_BE_SAME_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_BE_SAME_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Data type constraint.
// \ingroup constraints
//
// In case the two types \a A and \a B are not the same (ignoring all cv-qualifiers of both data
// types), a compilation error is created. The following example illustrates the behavior of this
// constraint:

   \code
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( double, double );        // No compilation error
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( double, const double );  // No compilation error (only cv-qualifiers differ)
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( double, float );         // Compilation error, different data types!
   \endcode

// In case the cv-qualifiers should not be ignored (e.g. 'double' and 'const double' should be
// considered to be unequal), use the blaze::BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE constraint.
*/
#define BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE(A,B) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_BE_SAME_TYPE_FAILED< ::blaze::IsSame<A,B>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_SAME_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_BE_SAME_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_BE_SAME_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_BE_SAME_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Data type constraint.
// \ingroup constraints
//
// In case the two types \a A and \a B are the same (ignoring all cv-qualifiers of both data
// types), a compilation error is created. The following example illustrates the behavior of
// this constraint:

   \code
   BLAZE_CONSTRAINT_MUST_NOT_BE_SAME_TYPE( double, float );         // No compilation error, different data types
   BLAZE_CONSTRAINT_MUST_NOT_BE_SAME_TYPE( double, const double );  // Compilation error (only cv-qualifiers differ)
   BLAZE_CONSTRAINT_MUST_NOT_BE_SAME_TYPE( double, double );        // Compilation error, same data type!
   \endcode

// In case the cv-qualifiers should not be ignored (e.g. 'double' and 'const double' should
// be considered to be unequal), use the blaze::BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE
// constraint.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_BE_SAME_TYPE(A,B) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_BE_SAME_TYPE_FAILED< !::blaze::IsSame<A,B>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_BE_SAME_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_BE_STRICTLY_SAME_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Data type constraint.
// \ingroup constraints
//
// In case the two types \a A and \a B are not the same, a compilation error is created. Note
// that this constraint even considers two types as unequal if the cv-qualifiers differ, e.g.

   \code
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( double, double );        // No compilation error
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( double, const double );  // Compilation error, different cv-qualifiers!
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( double, float );         // Compilation error, different data types!
   \endcode

// In case the cv-qualifiers should be ignored (e.g. 'double' and 'const double' should be
// considered to be equal), use the blaze::BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE constraint.
*/
#define BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE(A,B) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE_FAILED< ::blaze::IsStrictlySame<A,B>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_BE_STRICTLY_SAME_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Data type constraint.
// \ingroup constraints
//
// In case the two types \a A and \a B are the same, a compilation error is created. Note that
// this constraint even considers two types as unequal if the cv-qualifiers differ, e.g.

   \code
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE( double, float );         // No compilation error, different data types
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE( double, const double );  // No compilation error, different cv-qualifiers!
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE( double, double );        // Compilation error, same data type!
   \endcode

// In case the cv-qualifiers should be ignored (e.g. 'double' and 'const double' should be
// considered to be equal), use the blaze::BLAZE_CONSTRAINT_MUST_NOT_BE_SAME_TYPE constraint.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE(A,B) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE_FAILED< !::blaze::IsStrictlySame<A,B>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_BE_STRICTLY_SAME_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
