//=================================================================================================
/*!
//  \file blaze/math/simd/DerivedTypes.h
//  \brief Header file for the derived SIMD types
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

#ifndef _BLAZE_MATH_SIMD_DERIVEDTYPES_H_
#define _BLAZE_MATH_SIMD_DERIVEDTYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  DERIVED SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The SIMD data type for 'char'.
// \ingroup simd
*/
typedef SIMDTrait<char>::Type  simd_char_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'signed char'.
// \ingroup simd
*/
typedef SIMDTrait<signed char>::Type  simd_schar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned char'.
// \ingroup simd
*/
typedef SIMDTrait<unsigned char>::Type  simd_uchar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'wchar_t'.
// \ingroup simd
*/
typedef SIMDTrait<wchar_t>::Type  simd_wchar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<char>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<char> >::Type  simd_cchar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<signed char>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<signed char> >::Type  simd_cschar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned char>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<unsigned char> >::Type  simd_cuchar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<wchar_t>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<wchar_t> >::Type  simd_cwchar_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'short'.
// \ingroup simd
*/
typedef SIMDTrait<short>::Type  simd_short_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned short'.
// \ingroup simd
*/
typedef SIMDTrait<unsigned short>::Type  simd_ushort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<short>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<short> >::Type  simd_cshort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned short>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<unsigned short> >::Type  simd_cushort_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'int'.
// \ingroup simd
*/
typedef SIMDTrait<int>::Type  simd_int_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned int'.
// \ingroup simd
*/
typedef SIMDTrait<unsigned int>::Type  simd_uint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<int>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<int> >::Type  simd_cint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned int>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<unsigned int> >::Type  simd_cuint_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'long int'.
// \ingroup simd
*/
typedef SIMDTrait<long>::Type  simd_long_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned long int'.
// \ingroup simd
*/
typedef SIMDTrait<unsigned long>::Type  simd_ulong_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<long int>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<long> >::Type  simd_clong_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned long int>'.
// \ingroup simd
*/
typedef SIMDTrait< complex<unsigned long> >::Type  simd_culong_t;
//*************************************************************************************************

} // namespace blaze

#endif
