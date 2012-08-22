//=================================================================================================
/*!
//  \file blaze/util/Limits.h
//  \brief Numerical limits of built-in data types
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

#ifndef _BLAZE_UTIL_LIMITS_H_
#define _BLAZE_UTIL_LIMITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <limits>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Numerical limits of built-in data types.
// \ingroup util
//
// The Limits class provides numerical limits for the following built-in data types:
//
// <ul>
//    <li>Integral data types</li>
//    <ul>
//       <li>unsigned char, signed char, char, wchar_t</li>
//       <li>unsigned short, short</li>
//       <li>unsigned int, int</li>
//       <li>unsigned long, long</li>
//       <li>size_t, ptrdiff_t (for certain 64-bit compilers)</li>
//    </ul>
//    <li>Floating point data types</li>
//    <ul>
//       <li>float</li>
//       <li>double</li>
//       <li>long double</li>
//    </ul>
// </ul>
//
// Depending on the data type, the following limits can be used:
//
// - \b inf: The \a inf function is defined for all built-in data types. It returns the largest
//      possible positive value of the according data type.
// - \b ninf: The \a ninf function is defined for all signed integral and all floating point
//      data types. It returns the largest possible negative value of the according data type.
// - \b epsilon: The \a epsilon function is defined for all floating point data types and
//      returns the smallest possible difference between two values of the according data type.
// - \b accuracy: The \a accuracy function is defined for all floating point data types and
//      returns the computation accuracy of the corresponding data type. Due to the limited
//      floating point accuracy of a CPU this value is needed as computation threshold. This
//      value is used in most computations throughout the Blaze library.
// - \b fpuAccuracy: The \a fpuAccuracy function is defined for all floating point data types
//      and returns the floating point accuracy of the according point data type. Due to the
//      limited floating point accuracy of a CPU this value is needed as zero threshold in
//      computations.
//
// Code examples:

   \code
   // Positiv infinity value
   unsigned int ui = Limits<unsigned int>::inf();

   // Negative infinity value
   double d = Limits<double>::ninf();
   \endcode
*/
template< typename Type >
struct Limits
{};
//*************************************************************************************************




//=================================================================================================
//
//  SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<unsigned char> specialization.
// \ingroup util
*/
template<>
struct Limits<unsigned char>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive unsigned char value. */
   static inline unsigned char inf() { return std::numeric_limits<unsigned char>::max(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<char> specialization.
// \ingroup util
*/
template<>
struct Limits<char>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive char value. */
   static inline char inf () { return std::numeric_limits<char>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative char value. */
   static inline char ninf() { return std::numeric_limits<char>::min(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<signed char> specialization.
// \ingroup util
*/
template<>
struct Limits<signed char>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive signed char value. */
   static inline signed char inf () { return std::numeric_limits<signed char>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative signed char value. */
   static inline signed char ninf() { return std::numeric_limits<signed char>::min(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<wchar_t> specialization.
// \ingroup util
*/
template<>
struct Limits<wchar_t>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive wchar_t value. */
   static inline wchar_t inf () { return std::numeric_limits<wchar_t>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative wchar_t value. */
   static inline wchar_t ninf() { return std::numeric_limits<wchar_t>::min(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<unsigned short> specialization.
// \ingroup util
*/
template<>
struct Limits<unsigned short>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive unsigned short value. */
   static inline unsigned short inf() { return std::numeric_limits<unsigned short>::max(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<short> specialization.
// \ingroup util
*/
template<>
struct Limits<short>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive short value. */
   static inline short inf () { return std::numeric_limits<short>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative short value. */
   static inline short ninf() { return std::numeric_limits<short>::min(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<unsigned int> specialization.
// \ingroup util
*/
template<>
struct Limits<unsigned int>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive unsigned int value. */
   static inline unsigned int inf() { return std::numeric_limits<unsigned int>::max(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<int> specialization.
// \ingroup util
*/
template<>
struct Limits<int>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive int value. */
   static inline int inf () { return std::numeric_limits<int>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative int value. */
   static inline int ninf() { return std::numeric_limits<int>::min(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<unsigned long> specialization.
// \ingroup util
*/
template<>
struct Limits<unsigned long>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive unsigned long value. */
   static inline unsigned long inf() { return std::numeric_limits<unsigned long>::max(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<long> specialization.
// \ingroup util
*/
template<>
struct Limits<long>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive long value. */
   static inline long inf () { return std::numeric_limits<long>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative long value. */
   static inline long ninf() { return std::numeric_limits<long>::min(); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
#if defined(_WIN64)
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<size_t> specialization.
// \ingroup util
*/
template<>
struct Limits<std::size_t>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive size_t value. */
   static inline size_t inf() { return std::numeric_limits<size_t>::max(); }
};
/*! \endcond */
#endif
//*************************************************************************************************


//*************************************************************************************************
#if defined(_WIN64)
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<ptrdiff_t> specialization.
// \ingroup util
*/
template<>
struct Limits<ptrdiff_t>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive ptrdiff_t value. */
   static inline ptrdiff_t inf () { return std::numeric_limits<ptrdiff_t>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative ptrdiff_t value. */
   static inline ptrdiff_t ninf() { return std::numeric_limits<ptrdiff_t>::min(); }
};
/*! \endcond */
#endif
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<float> specialization.
// \ingroup util
*/
template<>
struct Limits<float>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive float value. */
   static inline float inf () { return  std::numeric_limits<float>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative float value. */
   static inline float ninf() { return -std::numeric_limits<float>::max(); }

   /*!\brief Machine epsilon.
   // \return The smallest possible difference between two float values. */
   static inline float epsilon() { return std::numeric_limits<float>::epsilon(); }

   /*!\brief The compuation accuracy of the Blaze library.
   // \return The computation threshold for single precision floating point values. */
   static inline float accuracy() { return 1E-6F; }

   /*!\brief The machine floating point accuracy.
   // \return The machine accuracy for single precision floating point values. */
   static inline float fpuAccuracy() { return 1E-12F; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<double> specialization.
// \ingroup util
*/
template<>
struct Limits<double>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive double value. */
   static inline double inf () { return  std::numeric_limits<double>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative double value. */
   static inline double ninf() { return -std::numeric_limits<double>::max(); }

   /*!\brief Machine epsilon.
   // \return The smallest possible difference between two double values. */
   static inline double epsilon() { return std::numeric_limits<double>::epsilon(); }

   /*!\brief The compuation accuracy of the Blaze library.
   // \return The computation threshold for double precision floating point values. */
   static inline double accuracy() { return 1E-8; }

   /*!\brief The machine floating point accuracy.
   // \return The machine accuracy for double precision floating point values. */
   static inline double fpuAccuracy() { return 1E-15; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Limits<long double> specialization.
// \ingroup util
*/
template<>
struct Limits<long double>
{
   /*!\brief Positive infinity value.
   // \return The largest possible positive long double value. */
   static inline long double inf () { return  std::numeric_limits<long double>::max(); }

   /*!\brief Negative infinity value.
   // \return The largest possible negative long double value. */
   static inline long double ninf() { return -std::numeric_limits<long double>::max(); }

   /*!\brief Machine epsilon.
   // \return The smallest possible difference between two long double values. */
   static inline long double epsilon() { return std::numeric_limits<long double>::epsilon(); }

   /*!\brief The compuation accuracy of the Blaze library.
   // \return The computation threshold for long double floating point values. */
   static inline long double accuracy() { return 1E-10L; }

   /*!\brief The machine floating point accuracy.
   // \return The machine accuracy for long double floating point values. */
   static inline long double fpuAccuracy() { return 1E-15L; }
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
