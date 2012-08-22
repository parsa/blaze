//=================================================================================================
/*!
//  \file blaze/system/Thresholds.h
//  \brief Header file for the thresholds for matrix/vector and matrix/matrix multiplications
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

#ifndef _BLAZE_SYSTEM_THRESHOLDS_H_
#define _BLAZE_SYSTEM_THRESHOLDS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/StaticAssert.h>




//=================================================================================================
//
//  THRESHOLDS
//
//=================================================================================================

#include <blaze/config/Thresholds.h>




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( blaze::DMATDVECMULT_THRESHOLD   > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDMATDVECMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDVECDMATMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDVECTDMATMULT_THRESHOLD > 0UL );
BLAZE_STATIC_ASSERT( blaze::DMATDMATMULT_THRESHOLD   > 0UL );
BLAZE_STATIC_ASSERT( blaze::DMATTDMATMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDMATDMATMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDMATTDMATMULT_THRESHOLD > 0UL );

}
/*! \endcond */
//*************************************************************************************************

#endif
