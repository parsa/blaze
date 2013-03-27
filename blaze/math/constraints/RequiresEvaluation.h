//=================================================================================================
/*!
//  \file blaze/math/constraints/RequiresEvaluation.h
//  \brief Constraint on the data type
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

#ifndef _BLAZE_MATH_CONSTRAINTS_REQUIRESEVALUATION_H_
#define _BLAZE_MATH_CONSTRAINTS_REQUIRESEVALUATION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/util/constraints/ConstraintTest.h>
#include <blaze/util/Suffix.h>


namespace blaze {

//=================================================================================================
//
//  MUST_REQUIRE_EVALUATION CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_REQUIRE_EVALUATION_FAILED;
template<> struct CONSTRAINT_MUST_REQUIRE_EVALUATION_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T does not require an intermediate evaluation within composite
// expressions, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_REQUIRE_EVALUATION(T) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_REQUIRE_EVALUATION_FAILED< blaze::RequiresEvaluation<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_REQUIRE_EVALUATION_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_REQUIRE_EVALUATION CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION_FAILED;
template<> struct CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T requires an intermediate evaluation within composite
// expressions, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION(T) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION_FAILED< !blaze::RequiresEvaluation<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
