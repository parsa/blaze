//=================================================================================================
/*!
//  \file blaze/math/shims/Equal.h
//  \brief Header file for the equal shim
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

#ifndef _BLAZE_MATH_SHIMS_EQUAL_H_
#define _BLAZE_MATH_SHIMS_EQUAL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/util/Complex.h>


namespace blaze {

//=================================================================================================
//
//  EQUAL SHIM
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for two single precision floating point values.
// \ingroup math_shims
//
// \param a The left-hand side single precision floating point value.
// \param b The right-hand side single precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two single precision floating point numbers. Due to the
// limited machine accuracy, a direct comparison of two floating point numbers should be avoided.
// This function offers the possibility to compare two floating-point values with a certain
// accuracy margin.
*/
inline bool equal( float a, float b )
{
   // Computing the absolute error
   if( std::fabs( a - b ) <= 1E-6F )
      return true;

   // Computing the relative error
   float relativeError;
   if( std::fabs(b) > std::fabs(a) )
      relativeError = std::fabs( ( a - b ) / b );
   else
      relativeError = std::fabs( ( a - b ) / a );

   if( relativeError <= 5E-4F )
      return true;
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for a single precision and a double precision floating point value.
// \ingroup math_shims
//
// \param a The left-hand side single precision floating point value.
// \param b The right-hand side double precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of a single precision and a double precision floating point
// number. Due to the limited machine accuracy, a direct comparison of two floating point numbers
// should be avoided. This function offers the possibility to compare two floating-point values
// with a certain accuracy margin.
*/
inline bool equal( float a, double b )
{
   return equal( a, static_cast<float>( b ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for a single precision and an extended precision floating point value.
// \ingroup math_shims
//
// \param a The left-hand side single precision floating point value.
// \param b The right-hand side extended precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of a single precision and an extended precision floating point
// number. Due to the limited machine accuracy, a direct comparison of two floating point numbers
// should be avoided. This function offers the possibility to compare two floating-point values
// with a certain accuracy margin.
*/
inline bool equal( float a, long double b )
{
   return equal( a, static_cast<float>( b ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for a double precision and a single precision floating point value.
// \ingroup math_shims
//
// \param a The left-hand side double precision floating point value.
// \param b The right-hand side single precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of a double precision and a single precision floating point
// number. Due to the limited machine accuracy, a direct comparison of two floating point numbers
// should be avoided. This function offers the possibility to compare two floating-point values
// with a certain accuracy margin.
*/
inline bool equal( double a, float b )
{
   return equal( static_cast<float>( a ), b );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for two double precision floating point values.
// \ingroup math_shims
//
// \param a The left-hand side double precision floating point value.
// \param b The right-hand side double precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two double precision floating point numbers. Due to the
// limited machine accuracy, a direct comparison of two floating point numbers should be avoided.
// This function offers the possibility to compare two floating-point values with a certain
// accuracy margin.
*/
inline bool equal( double a, double b )
{
   return std::fabs( a - b ) <= ( 1E-8 * std::fabs( a ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for a double precision and an extended precision floating point value.
// \ingroup math_shims
//
// \param a The left-hand side double precision floating point value.
// \param b The right-hand side extended precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of a double precision and an extended precision floating point
// number. Due to the limited machine accuracy, a direct comparison of two floating point numbers
// should be avoided. This function offers the possibility to compare two floating-point values
// with a certain accuracy margin.
*/
inline bool equal( double a, long double b )
{
   return std::fabs( a - b ) <= ( 1E-8L * std::fabs( b ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for an extended precision and a single precision floating point value.
// \ingroup math_shims
//
// \param a The left-hand side extended precision floating point value.
// \param b The right-hand side single precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of an extended precision and a single precision floating point
// number. Due to the limited machine accuracy, a direct comparison of two floating point numbers
// should be avoided. This function offers the possibility to compare two floating-point values
// with a certain accuracy margin.
*/
inline bool equal( long double a, float b )
{
   return equal( static_cast<float>( a ), b );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for an extended precision and a double precision floating point value.
// \ingroup math_shims
//
// \param a The left-hand side extended precision floating point value.
// \param b The right-hand side double precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of an extended precision and a double precision floating point
// number. Due to the limited machine accuracy, a direct comparison of two floating point numbers
// should be avoided. This function offers the possibility to compare two floating-point values
// with a certain accuracy margin.
*/
inline bool equal( long double a, double b )
{
   return std::fabs( a - b ) <= ( 1E-8L * std::fabs( a ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for two long double precision floating point values.
// \ingroup math_shims
//
// \param a The left-hand side extended precision floating point value.
// \param b The right-hand side extended precision floating point value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two long double precision floating point numbers. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should be
// avoided. This function offers the possibility to compare two floating-point values with a
// certain accuracy margin.
*/
inline bool equal( long double a, long double b )
{
   return std::fabs( a - b ) <= ( 1E-8L * std::fabs( a ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for a complex and a scalar value.
// \ingroup math_shims
//
// \param a The left-hand side complex value.
// \param b The right-hand side scalar value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of a complex and a scalar value. The function compares the
// real part of the complex value with the scalar. In case these two values match and in case
// the imaginary part is zero, the function returns \a true. Otherwise it returns \a false.
*/
template< typename T1    // Type of the left-hand side complex value
        , typename T2 >  // Type of the right-hand side scalar value
inline bool equal( complex<T1> a, T2 b )
{
   return equal( a.real(), b ) && isDefault( a.imag() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for a scalar and a complex value.
// \ingroup math_shims
//
// \param a The left-hand side scalar value.
// \param b The right-hand side complex value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of a scalar and a complex value. The function compares the
// scalar with the real part of the complex value. In case these two values match and in case
// the imaginary part is zero, the function returns \a true. Otherwise it returns \a false.
*/
template< typename T1    // Type of the left-hand side scalar value
        , typename T2 >  // Type of the right-hand side complex value
inline bool equal( T1 a, complex<T2> b )
{
   return equal( a, b.real() ) && isDefault( b.imag() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check for two complex values.
// \ingroup math_shims
//
// \param a The left-hand side complex value.
// \param b The right-hand side complex value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two complex numbers. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point values with a certain accuracy margin.
*/
template< typename T1    // Type of the left-hand side complex value
        , typename T2 >  // Type of the right-hand side complex value
inline bool equal( complex<T1> a, complex<T2> b )
{
   return equal( a.real(), b.real() ) && equal( a.imag(), b.imag() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generic equality check.
// \ingroup math_shims
//
// \param a First value/object.
// \param b Second value/object.
// \return \a true if the two values/objects are equal, \a false if not.
//
// The equal shim represents an abstract interface for testing two values/objects for equality.
// In case the two values/objects are equal, the function returns \a true, otherwise it returns
// \a false. Per default, the comparison of the two values/objects uses the equality operator
// operator==(). For built-in floating point data types a special comparison is selected that
// takes the limited machine accuracy into account.
*/
template< typename T1    // Type of the left-hand side value/object
        , typename T2 >  // Type of the right-hand side value/object
inline bool equal( const T1& a, const T2& b )
{
   return a == b;
}
//*************************************************************************************************

} // namespace blaze

#endif
