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
#include <blaze/util/Complex.h>
#include <blaze/util/NonCreatable.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename > class Rand;




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup random Random number generation
// \ingroup util
//
// The random number module provides the functionality to generate pseudo-random numbers within
// the Blaze library. In order to create series of random numbers, the following functions are
// are provided:
//
// - blaze::rand<T>();
// - blaze::getSeed();
// - blaze::setSeed( uint32_t seed );
//
// The templated rand() function is capable of generating random numbers for built-in integer
// and floating point data types as well as complex values. The rand() function can be given up
// to five parameters that depending on the type of the random number may have different meaning.
// The following example demonstrates the random number generation:

   \code
   using namespace blaze;

   // The setSeed function sets the seed for the random number generator. This function can
   // be used to set a specific seed for the random number generation, e.g. in order to
   // reproduce an exact series of random numbers. If the setSeed function is not used, the
   // random number generation uses a random seed.
   setSeed( 12345 );

   // In order to acquire the currently used seed, the getSeed() function can be used.
   uint32_t seed = getSeed();

   // Generating random numbers of built-in type. In case no range is provided, random numbers
   // of integral type are generated in the range [0..max], where max is the largest possible
   // value of the specified type and random floating point numbers are generated in the range
   // [0..1). In the example, the random double precision floating point value is created in
   // the range [2..4].
   int    i = rand<int>();
   double d = rand<double>( 2.0, 4.0 );

   // Generating random complex numbers. In case no range is provided, both the real and the
   // imaginary part are created depending on the element type of the complex number. In case
   // one range is specified, both the real and the imaginary part are created within the
   // specified range. In the example, both parts are created within the range [1..4]. If two
   // ranges are provided, the real part is created in the first range and the imaginary part
   // is created in the second range. The last example demonstrates this by restricting the
   // real part to the range [2..3] and the imaginary part to the range [1..5].
   complex<float> c1 = rand< complex<float> >();
   complex<float> c2 = rand< complex<float> >( 1.0F, 4.0F );
   complex<float> c3 = rand< complex<float> >( 2.0F, 3.0F, 1.0F, 5.0F );
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
   template< typename T > friend class Rand;
//    template< typename T > friend T        rand();
//    template< typename T > friend T        rand( T min, T max );
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
//  CLASS RAND
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default implementation of the Rand class for integral data types.
// \ingroup random
//
// This default implementation of the Rand class creates random, integral numbers in the range
// \f$ [0..max] \f$, where \a max is the maximal value of the given data type \a T.
*/
template< typename T >  // Type of the random number
class Rand
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   explicit inline Rand( T min, T max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator T() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   T value_;  //!< The random number.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor of the default implementation of the Rand class.
//
// The default constructor of the default implementation of the Rand class creates a random
// number in the range \f$ [0..max] \f$, where \a max is the maximal value of the given data
// type \a T.
*/
template< typename T >  // Type of the random number
inline Rand<T>::Rand()
   : value_()  // The random number
{
   boost::uniform_int<T> dist( 0, std::numeric_limits<T>::max() );
   value_ = dist( Random<RNG>::rng_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Range constructor of the default implementation of the Rand class.
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
//
// This constructor creates a random number in the range \f$ [min..max] \f$, where \a min must be
// smaller or equal to \a max. This requirement is only checked in debug mode. In release mode,
// no check is performed to enforce the validity of the values. Therefore the returned value is
// undefined if \a min is larger than \a max.
*/
template< typename T >  // Type of the random number
inline Rand<T>::Rand( T min, T max )
   : value_()  // The random number
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max value pair" );
   boost::uniform_smallint<T> dist( min, max );
   value_ = dist( Random<RNG>::rng_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion to the created random number.
//
// \return The random value.
*/
template< typename T >  // Type of the random number
inline Rand<T>::operator T() const
{
   return value_;
}
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION (FLOAT)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for single precision floating point values.
// \ingroup random
//
// This specialization of the Rand class creates random, single precision values in the range
// \f$ [0..1) \f$.
*/
template<>
class Rand<float>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   explicit inline Rand( float min, float max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator float() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   float value_;  //!< The random, single precision value.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor of the Rand<float> specialization.
//
// The default constructor of the Rand specialization creates a random number in the range
// \f$ [0..1) \f$.
*/
inline Rand<float>::Rand()
   : value_()  // The random, single precision value.
{
   boost::uniform_real<float> dist( 0.0, 1.0 );
   value_ = dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand<float> specialization.
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
//
// This constructor creates a random number in the range \f$ [min..max] \f$, where \a min must be
// smaller or equal to \a max. This requirement is only checked in debug mode. In release mode,
// no check is performed to enforce the validity of the values. Therefore the returned value is
// undefined if \a min is larger than \a max.
*/
inline Rand<float>::Rand( float min, float max )
   : value_()  // TODO
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max values" );
   boost::uniform_real<float> dist( min, max );
   value_ = dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random single precision number.
//
// \return The random value.
*/
inline Rand<float>::operator float() const
{
   return value_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION (DOUBLE)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for double precision floating point values.
// \ingroup random
//
// This specialization of the Rand class creates random, double precision values in the range
// \f$ [0..1) \f$.
*/
template<>
class Rand<double>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   explicit inline Rand( double min, double max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator double() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   double value_;  //!< The random, double precision value.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor of the Rand<double> specialization.
//
// The default constructor of the Rand specialization creates a random number in the range
// \f$ [0..1) \f$.
*/
inline Rand<double>::Rand()
   : value_()  // The random, double precision value.
{
   boost::uniform_real<double> dist( 0.0, 1.0 );
   value_ = dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand<double> specialization.
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
//
// This constructor creates a random number in the range \f$ [min..max] \f$, where \a min must be
// smaller or equal to \a max. This requirement is only checked in debug mode. In release mode,
// no check is performed to enforce the validity of the values. Therefore the returned value is
// undefined if \a min is larger than \a max.
*/
inline Rand<double>::Rand( double min, double max )
   : value_()  // TODO
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max values" );
   boost::uniform_real<double> dist( min, max );
   value_ = dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random double precision number.
//
// \return The random value.
*/
inline Rand<double>::operator double() const
{
   return value_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION (LONG DOUBLE)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for extended precision floating point values.
// \ingroup random
//
// This specialization of the Rand class creates random, extended precision values in the range
// \f$ [0..1) \f$.
*/
template<>
class Rand<long double>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   explicit inline Rand( long double min, long double max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator long double() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   long double value_;  //!< The random, extended precision value.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor of the Rand<long double> specialization.
//
// The default constructor of the Rand specialization creates a random number in the range
// \f$ [0..1) \f$.
*/
inline Rand<long double>::Rand()
   : value_()  // The random, extended precision value.
{
   boost::uniform_real<long double> dist( 0.0, 1.0 );
   value_ = dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand<long double> specialization.
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
//
// This constructor creates a random number in the range \f$ [min..max] \f$, where \a min must be
// smaller or equal to \a max. This requirement is only checked in debug mode. In release mode,
// no check is performed to enforce the validity of the values. Therefore the returned value is
// undefined if \a min is larger than \a max.
*/
inline Rand<long double>::Rand( long double min, long double max )
   : value_()  // TODO
{
   BLAZE_INTERNAL_ASSERT( min <= max, "Invalid min/max values" );
   boost::uniform_real<long double> dist( min, max );
   value_ = dist( Random<RNG>::rng_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random extended precision number.
//
// \return The random value.
*/
inline Rand<long double>::operator long double() const
{
   return value_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION (COMPLEX)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for complex values.
// \ingroup random
//
// This specialization of the Rand class creates random, complex values.
*/
template< typename T >  // Type of the values
class Rand< complex<T> >
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   explicit inline Rand( const T& min, const T& max );
   explicit inline Rand( const T& realmin, const T& realmax, const T& imagmin, const T& imagmax );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator complex<T>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   complex<T> value_;  //!< The random, complex value.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor of the Rand< complex<T> > specialization.
//
// The default constructor of the Rand specialization creates a random number in the range
// \f$ [0..1) \f$.
*/
template< typename T >  // Type of the values
inline Rand< complex<T> >::Rand()
   : value_( Rand<T>(), Rand<T>() )  // The random, complex value.
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand< complex<T> > specialization.
//
// \param min The smallest possible random value.
// \param max The largest possible random value.
//
// This constructor creates a random, complex number, where both the real and the imaginary part
// are in the range \f$ [min..max] \f$. Note that \a min must be smaller or equal to \a max. This
// requirement is only checked in debug mode. In release mode, no check is performed to enforce
// the validity of the values. Therefore the returned value is undefined if \a min is larger than
// \a max.
*/
template< typename T >  // Type of the values
inline Rand< complex<T> >::Rand( const T& min, const T& max )
   : value_( Rand<T>( min, max), Rand<T>( min, max ) )  // The random, complex value
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand< complex<T> > specialization.
//
// \param realmin The smallest possible random value for the real part.
// \param realmax The largest possible random value for the real part
// \param realmin The smallest possible random value for the imaginary part.
// \param realmax The largest possible random value for the imaginary part.
//
// This constructor creates a random, complex number, where the real part is in the range
// \f$ [realmin..realmax] \f$ and the imaginary part is in the range \f$ [imagmin..imagmax] \f$.
// \a realmin must be smaller or equal to \a realmax and \a imagmin must be smaller or equal to
// \a imagmax. These requirements are only checked in debug mode. In release mode, no check is
// performed to enforce the validity of the values. Therefore the returned value is undefined
// if \a realmin is larger than \a realmax or \a imagmin is larger than \a imagmax.
*/
template< typename T >  // Type of the values
inline Rand< complex<T> >::Rand( const T& realmin, const T& realmax, const T& imagmin, const T& imagmax )
   : value_( Rand<T>( realmin, realmax ), Rand<T>( imagmin, imagmax ) )  // The random, complex value
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random complex number.
//
// \return The random value.
*/
template< typename T >  // Type of the values
inline Rand< complex<T> >::operator complex<T>() const
{
   return value_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RANDOM NUMBER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Random number functions */
//@{
template< typename T >
inline T rand();

template< typename T, typename A1 >
inline T rand( const A1& a1 );

template< typename T, typename A1, typename A2 >
inline T rand( const A1& a1, const A2& a2 );

template< typename T, typename A1, typename A2, typename A3 >
inline T rand( const A1& a1, const A2& a2, const A3& a3 );

template< typename T, typename A1, typename A2, typename A3, typename A4 >
inline T rand( const A1& a1, const A2& a2, const A3& a3, const A4& a4 );

template< typename T, typename A1, typename A2, typename A3, typename A4, typename A5 >
inline T rand( const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5 );

inline uint32_t getSeed();
inline void     setSeed( uint32_t seed );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \return The generated random number.
//
// The rand() function returns a default random number depending on the given data type. In case
// of integral data types, the function returns a random number in the range \f$ [0..max] \f$,
// where \a max is the maximal value of the data type \a T. In case of floating point data types,
// the function returns a random number in the range \f$ [0..1) \f$. In case of complex data
// types, the function returns a random complex value where both the real and the imaginary part
// have been set according to the element type of the complex value (\f$ [0..max] \f$ for integral
// elements, \f$ [0..1) \f$ for floating point elements).
*/
template< typename T >  // Type of the random number
inline T rand()
{
   Rand<T> tmp;
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \param a1 First argument for the random number generation.
// \return The generated random number.
//
// This rand() function creates a random number based on the given argument \a a1.
*/
template< typename T     // Type of the random number
        , typename A1 >  // Type of the first argument
inline T rand( const A1& a1 )
{
   Rand<T> tmp( a1 );
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \param a1 First argument for the random number generation.
// \param a2 Second argument for the random number generation.
// \return The generated random number.
//
// This rand() function creates a random number based on the given arguments \a a1 and \a a2.
*/
template< typename T     // Type of the random number
        , typename A1    // Type of the first argument
        , typename A2 >  // Type of the second argument
inline T rand( const A1& a1, const A2& a2 )
{
   Rand<T> tmp( a1, a2 );
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \param a1 First argument for the random number generation.
// \param a2 Second argument for the random number generation.
// \param a3 Third argument for the random number generation.
// \return The generated random number.
//
// This rand() function creates a random number based on the given arguments \a a1, \a a2, and
// \a a3.
*/
template< typename T     // Type of the random number
        , typename A1    // Type of the first argument
        , typename A2    // Type of the second argument
        , typename A3 >  // Type of the third argument
inline T rand( const A1& a1, const A2& a2, const A3& a3 )
{
   Rand<T> tmp( a1, a2, a3 );
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \param a1 First argument for the random number generation.
// \param a2 Second argument for the random number generation.
// \param a3 Third argument for the random number generation.
// \param a4 Fourth argument for the random number generation.
// \return The generated random number.
//
// This rand() function creates a random number based on the given arguments \a a1, \a a2, \a a3,
// and \a a4.
*/
template< typename T     // Type of the random number
        , typename A1    // Type of the first argument
        , typename A2    // Type of the second argument
        , typename A3    // Type of the third argument
        , typename A4 >  // Type of the fourth argument
inline T rand( const A1& a1, const A2& a2, const A3& a3, const A4& a4 )
{
   Rand<T> tmp( a1, a2, a3, a4 );
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random number function.
// \ingroup random
//
// \param a1 First argument for the random number generation.
// \param a2 Second argument for the random number generation.
// \param a3 Third argument for the random number generation.
// \param a4 Fourth argument for the random number generation.
// \param a5 Fifth argument for the random number generation.
// \return The generated random number.
//
// This rand() function creates a random number based on the given arguments \a a1, \a a2, \a a3,
// \a a4, and \a a5.
*/
template< typename T     // Type of the random number
        , typename A1    // Type of the first argument
        , typename A2    // Type of the second argument
        , typename A3    // Type of the third argument
        , typename A4    // Type of the fourth argument
        , typename A5 >  // Type of the fifth argument
inline T rand( const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5 )
{
   Rand<T> tmp( a1, a2, a3, a4, a5 );
   return tmp;
}
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
