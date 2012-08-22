//=================================================================================================
/*!
//  \file blaze/util/valuetraits/IsOdd.h
//  \brief Header file for the IsEven value trait
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

#ifndef _BLAZE_UTIL_VALUETRAITS_ISODD_H_
#define _BLAZE_UTIL_VALUETRAITS_ISODD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check whether a compile time constant expression is odd.
// \ingroup value_traits
//
// This value trait tests whether the given integral value \a N is an odd value. In case the
// value is odd, the \a value member enumeration is set to 1, the nested type definition
// \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set
// to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsOdd<1>::value   // Evaluates to 1
   blaze::IsOdd<3>::Type    // Results in TrueType
   blaze::IsOdd<5>          // Is derived from TrueType
   blaze::IsOdd<2>::value   // Evaluates to 0
   blaze::IsOdd<4>::Type    // Results in FalseType
   blaze::IsOdd<6>          // Is derived from FalseType
   \endcode
*/
template< size_t N >
struct IsOdd : public SelectType<N%2,TrueType,FalseType>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = ( N%2 )?( 1 ):( 0 ) };
   typedef typename SelectType<N%2,TrueType,FalseType>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
