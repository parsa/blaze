//=================================================================================================
/*!
//  \file blaze/system/Restrict.h
//  \brief System settings for the restrict keyword
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

#ifndef _BLAZE_SYSTEM_RESTRICT_H_
#define _BLAZE_SYSTEM_RESTRICT_H_


//=================================================================================================
//
//  RESTRICT SETTINGS
//
//=================================================================================================

#include <blaze/config/Restrict.h>




//=================================================================================================
//
//  RESTRICT KEYWORD
//
//=================================================================================================

//*************************************************************************************************
/*!\def BLAZE_RESTRICT
// \brief Platform dependent setup of the restrict keyword.
// \ingroup system
*/
#if BLAZE_USE_RESTRICT

// Intel compiler
#  if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
#    define BLAZE_RESTRICT __restrict

// GNU compiler
#  elif defined(__GNUC__)
#    define BLAZE_RESTRICT __restrict

// Microsoft visual studio
#  elif defined(_MSC_VER)
#    define BLAZE_RESTRICT

// All other compilers
#  else
#    define BLAZE_RESTRICT

#  endif
#else
#  define BLAZE_RESTRICT
#endif
//*************************************************************************************************

#endif
