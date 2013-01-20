//=================================================================================================
/*!
//  \file blaze/system/Signature.h
//  \brief Header file for a compiler independent type/function signature macro.
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

#ifndef _BLAZE_SYSTEM_SIGNATURE_H_
#define _BLAZE_SYSTEM_SIGNATURE_H_


//=================================================================================================
//
//  SIGNATURE MACRO
//
//=================================================================================================

//*************************************************************************************************
/*!\def BLAZE_SIGNATURE
// \brief Platform dependent setup of the type/function signature macro.
// \ingroup system
//
// This macro contains the signature of the function the macro is used in. Note that the macro
// must only be used inside a function!
*/

// Intel compiler
#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
#  define BLAZE_SIGNATURE __PRETTY_FUNCTION__

// GNU compiler
#elif defined(__GNUC__)
#  define BLAZE_SIGNATURE __PRETTY_FUNCTION__

// Microsoft visual studio
#elif defined(_MSC_VER)
#  define BLAZE_SIGNATURE __FUNCSIG__

// All other compilers
#else
#  define BLAZE_SIGNATURE "Unknown function"
#endif
//*************************************************************************************************

#endif
