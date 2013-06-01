//=================================================================================================
/*!
//  \file blazetest/mathtest/RandomLimits.h
//  \brief Header file for the RandomLimits class template
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

#ifndef _BLAZETEST_MATHTEST_RANDOMLIMITS_H_
#define _BLAZETEST_MATHTEST_RANDOMLIMITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization limits for all built-in data types
//
// The RandomLimits class provides minimum and maximum limits for the random initialization of
// all numerical built-in data types. Via the \a min and \a max functions all numerical data
// types can be restricted to a fixed range of values. Integral data types can be restricted to
// the range \f$ [0..10] \f$, floating point data types to the range \f$ [0..1) \f$.
//
// Code examples:

   \code
   // Smallest allowed signed integer value
   int i = RandomLimits<int>::min();

   // Largest allowed double precision floating point value
   double d = RandomLimits<double>::max();
   \endcode
*/
template< typename Type >
struct RandomLimits
{
   /*!\brief Initialization minimum.
   // \return The smallest allowed initialization value. */
   static inline Type min() { return Type(0); }

   /*!\brief Initialization maximum.
   // \return The largest allowed initialization value. */
   static inline Type max() { return ( blaze::IsFloatingPoint<Type>::value )?( Type(1) ):( Type(10) ); }

   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( Type );
   /*! \endcond */
};
//*************************************************************************************************

} // namespace blazetest

#endif
