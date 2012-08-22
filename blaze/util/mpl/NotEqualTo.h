//=================================================================================================
/*!
//  \file blaze/util/mpl/NotEqualTo.h
//  \brief Header file for the NotEqualTo class template
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

#ifndef _BLAZE_UTIL_MPL_NOTEQUALTO_H_
#define _BLAZE_UTIL_MPL_NOTEQUALTO_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type value comparison.
// \ingroup mpl
//
// The NotEqualTo class templates evaluates whether the two given conditions \a T1 and \a T2
// are not equal to each other. In case \a T1::value is not equal to \a T2::value, the \a value
// member enumeration is set to 1, otherwise it is set to 0.

   \code
   EqualTo< IsDouble<double>, IsFloat<double> >::value  // Evaluates to 1
   EqualTo< IsDouble<float> , IsFloat<float>  >::value  // Evaluates to 1
   EqualTo< IsDouble<double>, IsFloat<float>  >::value  // Evaluates to 0
   EqualTo< IsDouble<float> , IsFloat<double> >::value  // Evaluates to 0
   \endcode
*/
template< typename T1    // Type of the left-hand side condition
        , typename T2 >  // Type of the right-hand side condition
struct NotEqualTo
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = ( T1::value && !T2::value ) || ( !T1::value && T2::value ) };
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
