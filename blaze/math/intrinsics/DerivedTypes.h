//=================================================================================================
/*!
//  \file blaze/math/intrinsics/DerivedTypes.h
//  \brief Header file for the derived intrinsic types
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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
typedef IntrinsicTrait<short>::Type  simd_short_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'unsigned short'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<unsigned short>::Type  simd_ushort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<int>::Type  simd_int_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'unsigned int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<unsigned int>::Type  simd_uint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'long int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<long>::Type  simd_long_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'unsigned long int'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait<unsigned long>::Type  simd_ulong_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'complex<short>'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait< complex<short> >::Type  simd_cshort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'complex<unsigned short>'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait< complex<unsigned short> >::Type  simd_cushort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'complex<int>'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait< complex<int> >::Type  simd_cint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'complex<unsigned int>'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait< complex<unsigned int> >::Type  simd_cuint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'complex<long int>'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait< complex<long> >::Type  simd_clong_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The intrinsic data type for 'complex<unsigned long int>'.
// \ingroup intrinsics
*/
typedef IntrinsicTrait< complex<unsigned long> >::Type  simd_culong_t;
//*************************************************************************************************

} // namespace blaze

#endif
