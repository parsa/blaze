//=================================================================================================
/*!
//  \file blaze/util/Unused.h
//  \brief Header file for the UNUSED_PARAMETER function template
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

#ifndef _BLAZE_UTIL_UNUSED_H_
#define _BLAZE_UTIL_UNUSED_H_


namespace blaze {

//=================================================================================================
//
//  UNUSED_PARAMETER FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Suppression of unused parameter warnings.
// \ingroup util
//
// \return void
//
// The UNUSED_PARAMETER function provides the functionality to suppress warnings about up to
// six unused parameters. Usually this problem occurs in case a parameter is given a name but
// is not used within the function:

   \code
   void f( int x )
   {}  // x is not used within f. This may result in an unused parameter warning.
   \endcode

// A possible solution is to keep the parameter unnamed:

   \code
   void f( int )
   {}  // No warning about unused parameter is issued
   \endcode

// However, there are situations where is approach is not possible, as for instance in case the
// variable must be documented via Doxygen. For these cases, the UNUSED_PARAMETER class can be
// used to suppress the warnings:

   \code
   void f( int x )
   {
      UNUSED_PARAMETER( x );  // Suppresses the unused parameter warnings
   }
   \endcode

// The NullType class represents an invalid or terminating data type for generic codes. For
// instance, the TypeList class uses the NullType as terminating data type for the type list.
*/
template< typename T1 >
inline void UNUSED_PARAMETER( const T1& )
{}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UNUSED_PARAMETER function for two parameters.
// \ingroup util
//
// \return void
*/
template< typename T1, typename T2 >
inline void UNUSED_PARAMETER( const T1&, const T2& )
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UNUSED_PARAMETER function for three parameters.
// \ingroup util
//
// \return void
*/
template< typename T1, typename T2, typename T3 >
inline void UNUSED_PARAMETER( const T1&, const T2&, const T3& )
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UNUSED_PARAMETER function for four parameters.
// \ingroup util
//
// \return void
*/
template< typename T1, typename T2, typename T3, typename T4 >
inline void UNUSED_PARAMETER( const T1&, const T2&, const T3&, const T4& )
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UNUSED_PARAMETER function for five parameters.
// \ingroup util
//
// \return void
*/
template< typename T1, typename T2, typename T3, typename T4, typename T5 >
inline void UNUSED_PARAMETER( const T1&, const T2&, const T3&, const T4&, const T5& )
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UNUSED_PARAMETER function for six parameters.
// \ingroup util
//
// \return void
*/
template< typename T1, typename T2, typename T3, typename T4, typename T5, typename T6 >
inline void UNUSED_PARAMETER( const T1&, const T2&, const T3&, const T4&, const T5&, const T6& )
{}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
