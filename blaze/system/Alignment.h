//=================================================================================================
/*!
//  \file blaze/system/Alignment.h
//  \brief System specific memory alignment definitions
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

#ifndef _BLAZE_SYSTEM_ALIGNMENT_H_
#define _BLAZE_SYSTEM_ALIGNMENT_H_


//=================================================================================================
//
//  ALIGNMENT MACRO
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\def BLAZE_ALIGN(n)
// \brief Platform dependent macro for the alignment of memory.
// \ingroup system
*/

// GNU, Intel, PGI, and IBM C++ compiler
#if (defined __GNUC__) || (defined __PGI) || (defined __IBMCPP__)
#  define BLAZE_ALIGN( n ) __attribute__((aligned(n)))

// Microsoft Visual C++ compiler
#elif (defined _MSC_VER)
#  define BLAZE_ALIGN( n ) __declspec(align(n))

// Default case for other compilers
#else
#  error Compiler-specific alignment undefined!
#endif

/*! \endcond */
//*************************************************************************************************

#endif
