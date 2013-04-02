//=================================================================================================
/*!
//  \file blaze/config/Vectorization.h
//  \brief Configuration of the vectorization policy of the Blaze library
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


//*************************************************************************************************
/*!\brief Compilation switch for (de-)activation of the Blaze vectorization.
// \ingroup config
//
// This compilation switch enables/disables vectorization of mathematical expressions via
// the SSE, AVX, and/or MIC instruction sets. In case the switch is set to 1 (i.e. in case
// vectorization is enabled) ,the Blaze library attempts to vectorize the linear algebra
// operations by SSE, AVX, and/or MIC intrinsics (depending on which instruction set is
// available on the target platform). In case the switch is set to 0 (i.e. vectorization
// is disabled), the Blaze library chooses default, non-vectorized functionality for the
// operations. Note that deactivating the vectorization may pose a severe performance
// limitation for a large number of operations!
//
// Possible settings for the vectorization switch:
//  - Deactivated: \b 0
//  - Activated  : \b 1 (default)
*/
#define BLAZE_USE_VECTORIZATION 1
//*************************************************************************************************
