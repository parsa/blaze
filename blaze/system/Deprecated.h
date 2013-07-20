//=================================================================================================
/*!
//  \file blaze/system/Deprecated.h
//  \brief System specific deprecated tags
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

#ifndef _BLAZE_SYSTEM_DEPRECATED_H_
#define _BLAZE_SYSTEM_DEPRECATED_H_


//=================================================================================================
//
//  DEPRECATED MACRO
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\def BLAZE_DEPRECATED(func)
// \brief Platform dependent macro for marking a function as deprecated.
// \ingroup system
*/

// GNU, Intel, PGI, and IBM C++ compiler
#if (defined __GNUC__) || (defined __PGI) || (defined __IBMCPP__)
#  define BLAZE_DEPRECATED( func ) func __attribute__((deprecated))

// Microsoft Visual C++ compiler
#elif (defined _MSC_VER)
#  define BLAZE_DEPRECATED( func ) __declspec(deprecated) func

// Default case for other compilers
#else
#  error Compiler-specific deprecated tag undefined!
#endif

/*! \endcond */
//*************************************************************************************************

#endif
