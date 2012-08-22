//=================================================================================================
/*!
//  \file blaze/system/BLAS.h
//  \brief System settings for the BLAS mode
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

#ifndef _BLAZE_SYSTEM_BLAS_H_
#define _BLAZE_SYSTEM_BLAS_H_


//=================================================================================================
//
//  BLAS MODE CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compilation switch for the BLAS mode.
// \ingroup system
//
// This compilation switch enables/disables the BLAS mode. In case the BLAS mode is enabled,
// several basic linear algebra functions (such as for instance matrix-matrix multiplications
// between two dense matrices) are handled by performance optimized BLAS functions. Note that
// in this case it is mandatory to provide the according BLAS header file for the compilation
// of the Blaze library. In case the BLAS mode is disabled, all linear algebra functions use
// the default implementations of the Blaze library and therefore BLAS is not a requirement
// for the compilation process.
//
// Possible settings for the BLAS switch:
//  - Deactivated: \b 0
//  - Activated  : \b 1
//
// Note that changing the setting of the BLAS mode requires a recompilation of the Blaze
// library. Also note that this switch is automatically set by the configuration script of
// the Blaze library.
*/
#define BLAZE_BLAS_MODE 0
//*************************************************************************************************

#endif
