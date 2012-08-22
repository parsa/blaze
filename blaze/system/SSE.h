//=================================================================================================
/*!
//  \file blaze/system/SSE.h
//  \brief System settings for the SSE mode
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

#ifndef _BLAZE_SYSTEM_SSE_H_
#define _BLAZE_SYSTEM_SSE_H_


//=================================================================================================
//
//  SSE MODE CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compilation switch for the SSE mode.
// \ingroup system
//
// This compilation switch enables/disables the SSE mode. In case the SSE mode is enabled
// (i.e. in case SSE functionality is available) the Blaze library attempts to vectorize
// the linear algebra operations by SSE intrinsics. In case the SSE mode is disabled, the
// Blaze library chooses default, non-vectorized functionality for the operations.
*/
#if defined(__SSE__) || ( _M_IX86_FP > 0 )
#  define BLAZE_SSE_MODE 1
#else
#  define BLAZE_SSE_MODE 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the SSE2 mode.
// \ingroup system
//
// This compilation switch enables/disables the SSE2 mode. In case the SSE2 mode is enabled
// (i.e. in case SSE2 functionality is available) the Blaze library attempts to vectorize
// the linear algebra operations by SSE2 intrinsics. In case the SSE2 mode is disabled, the
// Blaze library chooses default, non-vectorized functionality for the operations.
*/
#if defined(__SSE2__) || ( _M_IX86_FP > 1 )
#  define BLAZE_SSE2_MODE 1
#else
#  define BLAZE_SSE2_MODE 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the SSE3 mode.
// \ingroup system
//
// This compilation switch enables/disables the SSE3 mode. In case the SSE3 mode is enabled
// (i.e. in case SSE3 functionality is available) the Blaze library attempts to vectorize
// the linear algebra operations by SSE3 intrinsics. In case the SSE3 mode is disabled, the
// Blaze library chooses default, non-vectorized functionality for the operations.
*/
#if defined(__SSE3__)
#  define BLAZE_SSE3_MODE 1
#else
#  define BLAZE_SSE3_MODE 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the SSSE3 mode.
// \ingroup system
//
// This compilation switch enables/disables the SSSE3 mode. In case the SSSE3 mode is enabled
// (i.e. in case SSSE3 functionality is available) the Blaze library attempts to vectorize
// the linear algebra operations by SSSE3 intrinsics. In case the SSSE3 mode is disabled, the
// Blaze library chooses default, non-vectorized functionality for the operations.
*/
#if defined(__SSSE3__)
#  define BLAZE_SSSE3_MODE 1
#else
#  define BLAZE_SSSE3_MODE 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the SSE4 mode.
// \ingroup system
//
// This compilation switch enables/disables the SSE4 mode. In case the SSE4 mode is enabled
// (i.e. in case SSE4 functionality is available) the Blaze library attempts to vectorize
// the linear algebra operations by SSE4 intrinsics. In case the SSE4 mode is disabled,
// the Blaze library chooses default, non-vectorized functionality for the operations.
*/
#if defined(__SSE4_2__) || defined(__SSE4_1__)
#  define BLAZE_SSE4_MODE 1
#else
#  define BLAZE_SSE4_MODE 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the AVX mode.
// \ingroup system
//
// This compilation switch enables/disables the AVX mode. In case the AVX mode is enabled
// (i.e. in case AVX functionality is available) the Blaze library attempts to vectorize
// the linear algebra operations by AVX intrinsics. In case the AVX mode is disabled,
// the Blaze library chooses default, non-vectorized functionality for the operations.
*/
#if defined(__AVX__)
#  define BLAZE_AVX_MODE 1
#else
#  define BLAZE_AVX_MODE 0
#endif
//*************************************************************************************************




//=================================================================================================
//
//  SSE INCLUDE FILE CONFIGURATION
//
//=================================================================================================

#if BLAZE_AVX_MODE
#  include <immintrin.h>
#elif BLAZE_SSE4_MODE
#  include <smmintrin.h>
#elif BLAZE_SSSE3_MODE
#  include <tmmintrin.h>
#elif BLAZE_SSE3_MODE
#  include <pmmintrin.h>
#elif BLAZE_SSE2_MODE
#  include <emmintrin.h>
#elif BLAZE_SSE_MODE
#  include <xmmintrin.h>
#endif

#endif
