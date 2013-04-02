//=================================================================================================
/*!
//  \file blaze/math/intrinsics/DerivedTypes.h
//  \brief Header file for the derived intrinsic types
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

#ifndef _BLAZE_MATH_INTRINSICS_DERIVEDTYPES_H_
#define _BLAZE_MATH_INTRINSICS_DERIVEDTYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/IntrinsicTrait.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  DERIVED INTRINSIC TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The intrinsic data type for 'short'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<short>::Type  sse_short_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'unsigned short'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<unsigned short>::Type  sse_ushort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<int>::Type  sse_int_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'unsigned int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<unsigned int>::Type  sse_uint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'long int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<long>::Type  sse_long_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'unsigned long int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<unsigned long>::Type  sse_ulong_t;
//*************************************************************************************************

} // namespace blaze

#endif
