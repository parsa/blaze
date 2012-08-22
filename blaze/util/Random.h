//=================================================================================================
/*!
//  \file blaze/util/Random.h
//  \brief Implementation of a random number generator.
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

#ifndef _BLAZE_UTIL_RANDOM_H_
#define _BLAZE_UTIL_RANDOM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <ctime>
#include <limits>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <blaze/system/Random.h>
#include <blaze/util/Assert.h>
#include <blaze/util/NonCreatable.h>
#include <blaze/util/Types.h>


namespace blaze {

//*************************************************************************************************
/*!\defgroup random Random number generation
// \ingroup util
//
// The random number module provides the functionality to generate pseudo-random numbers within
// the Blaze library. In order to create series of random numbers, the following four functions
// are provided:
//
// - blaze::rand<T>();
// - blaze::rand<T>( T min, T max );
// - blaze::getSeed();
// - blaze::setSeed( uint32_t seed );
//
// The templated rand() functions are capable of generate random numbers for both built-in integer
// and floating point data types. The following example demonstrates the random number generation:

   \code
   // The setSeed function sets the seed for the random number generator. This function can
   // be used to set a specific seed for the random number generation, e.g. in order to
   // reproduce an exact series of random numbers. If the setSeed function is not used, the
   // random number generation uses a random seed.
   setSeed( 12345 );

   // Generating random numbers
   int    i = rand<int>();
   float  f = rand<float>();
   double d = rand<double>();
   \endcode

// \b Note: In order to reproduce certain series of random numbers, the seed of the random number
// generator has to be set explicitly via the setSeed() function. Otherwise a random seed is used
// for the random number generation.
*/
/*!\brief Random number generator.
// \ingroup random
//
// The Random class encapsulates the initialization of the given random number generator with
// a pseudo-random seed obtained by the std::time() function. Currently, the mersenne-twister
// mt19937 as provided by the boost library is used per default. For more information see the
// class description of the boost library:
//
//   http://www.boost.org/doc/libs/1_35_0/libs/random/random-generators.html#mersenne_twister\n
//   http://www.boost.org/doc/libs/1_35_0/boost/random/mersenne_twister.hpp
*/
template< typename Type >  // Type of the random number generator
class Random : private NonCreatable
{
 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static uint32_t seed_;  //!< The current seed for the variate generator.
   static Type     rng_;   //!< The mersenne twister variate generator.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T > friend T        rand();
   template< typename T > friend T        rand( T min, T max );
                          friend uint32_t getSeed();
                          friend void     setSeed( uint32_t seed );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename Type > uint32_t Random<Type>::seed_( static_cast<uint32_t>( std::time(0) ) );
template< typename Type > Type     Random<Type>::rng_ ( seed_ );




//=================================================================================================
//
//  RANDOM NUMBER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Random number functions */
//@{
template< typename T > inline T    rand();
template< typename T > inline T    rand( T min, T max );
                       inline void setSeed( uint32_t seed );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \return The generated random number.
//
// The rand() function returns a random number depending on the given data type. In case of
// integral data types, the function returns a random number in the range \f$ [0..max] \f$,
// where \a max is the maximal value of the data type \a T. In case of floating point data
// types, the function returns a random number in the range \f$ [0..1) \f$.
*/
template< typename T >
inline T rand()
{
   boost::uniform_int<T> dist( 0, std::numeric_limits<T>::max() );
   return dist( Random<RNG>::rng_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the rand function for random float values.
// \ingroup random
//
// \return The generated random float value.
//
// The function returns random float values in the range \f$ [0..1) \f$.
*/
template<>
inline float rand<float>()
{
   boost::uniform_real<float> dist( 0.0F, 1.0F );
   return dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the rand function for random double values.
// \ingroup random
//
// \return The generated random double value.
//
// The function returns random double values in the range \f$ [0..1) \f$.
*/
template<>
inline double rand<double>()
{
   boost::uniform_real<double> dist( 0.0, 1.0 );
   return dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the rand function for random long double values.
// \ingroup random
//
// \return The generated random long double value.
//
// The function returns random long double values in the range \f$ [0..1) \f$.
*/
template<>
inline long double rand<long double>()
{
   boost::uniform_real<long double> dist( 0.0L, 1.0L );
   return dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
// \return The generated random number.
//
// This function returns a random number in the range \f$ [min..max] \f$, where \a min must be
// smaller or equal to \a max. This requirement is only checked in debug mode. In release mode,
// no check is performed to enforce the validity of the values. Therefore the returned value is
// undefined if \a min is larger than \a max.
*/
template< typename T >
inline T rand( T min, T max )
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max value pair" );
   boost::uniform_smallint<T> dist( min, max );
   return dist( Random<RNG>::rng_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the min/max-rand function for random float values.
// \ingroup random
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
// \return The generated random float value.
//
// The function returns random float values in the range \f$ [min..max) \f$, where \a min
// must be smaller or equal to \a max. This requirement is only checked in debug mode. In
// release mode, no check is performed to enforce the validity of the values. Therefore the
// returned value is undefined if \a min is larger than \a max.
*/
template<>
inline float rand<float>( float min, float max )
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max value pair" );
   boost::uniform_real<float> dist( min, max );
   return dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the min/max-rand function for random double values.
// \ingroup random
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
// \return The generated random double value.
//
// The function returns random double values in the range \f$ [min..max) \f$, where \a min
// must be smaller or equal to \a max. This requirement is only checked in debug mode. In
// release mode, no check is performed to enforce the validity of the values. Therefore the
// returned value is undefined if \a min is larger than \a max.
*/
template<>
inline double rand<double>( double min, double max )
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max values" );
   boost::uniform_real<double> dist( min, max );
   return dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the min/max-rand function for random long double values.
// \ingroup random
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
// \return The generated random long double value.
//
// The function returns random long double values in the range \f$ [min..max) \f$, where
// \a min must be smaller or equal to \a max. This requirement is only checked in debug
// mode. In release mode, no check is performed to enforce the validity of the values.
// Therefore the returned value is undefined if \a min is larger than \a max.
*/
template<>
inline long double rand<long double>( long double min, long double max )
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max values" );
   boost::uniform_real<long double> dist( min, max );
   return dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current seed of the random number generator.
// \ingroup random
//
// \return The current seed of the random number generator.
*/
inline uint32_t getSeed()
{
   return Random<RNG>::seed_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the seed of the random number generator.
// \ingroup random
//
// \param seed The new seed for the random number generator.
// \return void
//
// This function can be used to set the seed for the random number generation in order to
// create a reproducible series of random numbers.
*/
inline void setSeed( uint32_t seed )
{
   Random<RNG>::seed_ = seed;
   Random<RNG>::rng_.seed( seed );
}
//*************************************************************************************************

} // namespace blaze

#endif
