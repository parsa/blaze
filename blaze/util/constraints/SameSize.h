//=================================================================================================
/*!
//  \file blaze/util/constraints/SameSize.h
//  \brief Constraint on the size of two data types
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

#ifndef _BLAZE_UTIL_CONSTRAINTS_SAMESIZE_H_
#define _BLAZE_UTIL_CONSTRAINTS_SAMESIZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/ConstraintTest.h>
#include <blaze/util/Suffix.h>
#include <blaze/util/typetraits/HaveSameSize.h>


namespace blaze {

//=================================================================================================
//
//  MUST_HAVE_SAME_SIZE CONSTRAINT
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
template< bool > struct CONSTRAINT_MUST_HAVE_SAME_SIZE_FAILED;
template<> struct CONSTRAINT_MUST_HAVE_SAME_SIZE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of two data types.
// \ingroup constraints
//
// In case the types \a T1 and \a T2 don't have the same size, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE(T1,T2) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_HAVE_SAME_SIZE_FAILED< ::blaze::HaveSameSize<T1,T2>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_HAVE_SAME_SIZE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_HAVE_SAME_SIZE CONSTRAINT
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
template< bool > struct CONSTRAINT_MUST_NOT_HAVE_SAME_SIZE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_HAVE_SAME_SIZE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of two data types.
// \ingroup constraints
//
// In case the types \a T1 and \a T2 have the same size, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_HAVE_SAME_SIZE(T1,T2) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_HAVE_SAME_SIZE_FAILED< !::blaze::HaveSameSize<T1,T2>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_HAVE_SAME_SIZE_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
