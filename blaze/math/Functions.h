//=================================================================================================
/*!
//  \file blaze/math/Functions.h
//  \brief Header file for mathematical functions
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

#ifndef _BLAZE_MATH_FUNCTIONS_H_
#define _BLAZE_MATH_FUNCTIONS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <limits>
#include <blaze/math/MathTrait.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Integral.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  MATHEMATICAL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Mathematical utility functions */
//@{
template< typename T >
inline const T sign( T a );

template< typename T >
inline size_t digits( T a );

template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::HighType
   min( const T1& a, const T2& b );

template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::HighType, T3 >::HighType
   min( const T1& a, const T2& b, const T3& c );

template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::HighType
   max( const T1& a, const T2& b );

template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::HighType, T3 >::HighType
   max( const T1& a, const T2& b, const T3& c );

template< typename T >
inline T round( T a );

template< typename T1, typename T2 >
inline bool lessThan( T1 a, T2 b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sign function.
// \ingroup math
//
// \param a The signed value.
// \return 1 if the value is greater than or equal to zero, -1 if the value is smaller than zero.
//
// The sign function only works for signed built-in data types. The attempt to use unsigned data
// types or user-defined class types will result in a compile time error.
*/
template< typename T >
inline const T sign( T a )
{
   BLAZE_STATIC_ASSERT( std::numeric_limits<T>::is_signed );
   return ( a < T(0) )?( T(-1) ):( T(1) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of valid digits of an integral value.
// \ingroup math
//
// \param a The integral value.
// \return The number of valid digits.
//
// This function counts the number of valid digits in the given integral value.

   \code
   digits( 100   );  // Returns 3
   digits( 12345 );  // Returns 5
   digits( 0     );  // Returns 0
   \endcode

// The digits function only works for integral built-in data types. The attempt to use any
// other type will result in a compile time error.
*/
template< typename T >
inline size_t digits( T a )
{
   BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( T );

   size_t count( 0 );

   while( a != 0 ) {
      a /= 10;
      ++count;
   }

   return count;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Minimum function for two arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return The minimum of the two values.
//
// This function returns the minimum of the two given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
*/
template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::HighType min( const T1& a, const T2& b )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1 and T2 are equal, a reference return type could result in a
   //    bug if combined with literals
   //  - in case the two data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a < b )?( a ):( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Minimum function for three arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \param c Third value
// \return The minimum of the three values.
//
// This function returns the minimum of the three given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
*/
template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::HighType, T3 >::HighType
   min( const T1& a, const T2& b, const T3& c )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1, T2, and T3 are equal, a reference return type could result in
   //    a bug if combined with literals
   //  - in case the three data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a < b )?( ( a < c )?( a ):( c ) ):( ( b < c )?( b ):( c ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Maximum function for two arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return The maximum of the two values.
//
// This function returns the maximum of the two given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
*/
template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::HighType max( const T1& a, const T2& b )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1 and T2 are equal, a reference return type could result in a
   //    bug if combined with literals
   //  - in case the two data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a < b )?( b ):( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Maximum function for three arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \param c Third value.
// \return The maximum of the three values.
//
// This function returns the maximum of the three given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
*/
template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::HighType, T3 >::HighType
   max( const T1& a, const T2& b, const T3& c )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1, T2, and T3 are equal, a reference return type could result in
   //    a bug if combined with literals
   //  - in case the three data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a < b )?( ( b < c )?( c ):( b ) ):( ( a < c )?( c ):( a ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rounds the given input value.
// \ingroup math
//
// \param a Value to be rounded.
// \return The rounded value.
//
// This function rounds the given input value. In case the first digit after the comma
// is smaller than five the value is rounded down. Otherwise it is rounded up. Note that
// this function only works for integral and floating point types. The attempt to use the
// function for any other type will result in a compile time error.
*/
template< typename T >
inline T round( T a )
{
   BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( T );
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Rounds the given single precision floating point value.
// \ingroup math
//
// \param a Value to be rounded.
// \return The rounded value.
//
// This function rounds the given single precision floating point value. In case the first
// digit after the comma is smaller than five the value is rounded down. Otherwise it is
// rounded up.
*/
template<>
inline float round( float a )
{
   return std::floor( a + 0.5F );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Rounds the given double precision floating point value.
// \ingroup math
//
// \param a Value to be rounded.
// \return The rounded value.
//
// This function rounds the given double precision floating point value. In case the first
// digit after the comma is smaller than five the value is rounded down. Otherwise it is
// rounded up.
*/
template<>
inline double round<double>( double a )
{
   return std::floor( a + 0.5 );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Rounds the given long double precision floating point value.
// \ingroup math
//
// \param a Value to be rounded.
// \return The rounded value.
//
// This function rounds the given long double precision floating point value. In case the
// first digit after the comma is smaller than five the value is rounded down. Otherwise
// it is rounded up.
*/
template<>
inline long double round<long double>( long double a )
{
   return std::floor( a + 0.5L );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Less-than comparison for integral data types.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the first value is smaller than the second, \a false if not.
//
// Less-than function for the comparison of two integral values.
*/
template< typename T >
inline bool lessThan_backend( T a, T b )
{
   return a < b;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Less-than comparison for single precision floating point values.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the first value is smaller than the second, \a false if not.
//
// Less-than function for the comparison of two single precision floating point numbers. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should
// be avoided. This functions offers the possibility to compare two floating-point values with
// a certain accuracy margin.
*/
template<>
inline bool lessThan_backend<float>( float a, float b )
{
   return ( b - a ) > 1E-8F;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Less-than comparison for double precision floating point values.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the first value is smaller than the second, \a false if not.
//
// Less-than function for the comparison of two double precision floating point numbers. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should
// be avoided. This functions offers the possibility to compare two floating-point values with
// a certain accuracy margin.
*/
template<>
inline bool lessThan_backend<double>( double a, double b )
{
   return ( b - a ) > 1E-8;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Less-than comparison for long double precision floating point values.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the first value is smaller than the second, \a false if not.
//
// Less-than function for the comparison of two long double precision floating point numbers. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should be
// avoided. This functions offers the possibility to compare two floating-point values with a
// certain accuracy margin.
*/
template<>
inline bool lessThan_backend<long double>( long double a, long double b )
{
   return ( b - a ) > 1E-10;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generic less-than comparison.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the first value is smaller than the second, \a false if not.
//
// Generic less-than comparison between to numeric values. Depending on the types of the
// two arguments, a special comparison for floating point values is selected that takes
// the limited machine accuracy into account.
*/
template< typename T1, typename T2 >
inline bool lessThan( T1 a, T2 b )
{
   typedef typename MathTrait<T1,T2>::HighType  High;
   return lessThan_backend<High>( a, b );
}
//*************************************************************************************************

} // namespace blaze

#endif
