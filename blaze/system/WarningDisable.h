//=================================================================================================
/*!
//  \file blaze/system/WarningDisable.h
//  \brief Deactivation of compiler specific warnings
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

#ifndef _BLAZE_SYSTEM_WARNINGDISABLE_H_
#define _BLAZE_SYSTEM_WARNINGDISABLE_H_


//=================================================================================================
//
//  MICROSOFT VISUAL STUDIO WARNINGS
//
//=================================================================================================

#if defined(_MSC_VER) && (_MSC_VER >= 1400)

   // Disables a 'deprecated' warning for some standard library functions. This warning
   // is emitted when you use some perfectly conforming library functions in a perfectly
   // correct way, and also by some of Microsoft's own standard library code. For more
   // information about this particular warning, see
   // http://msdn.microsoft.com/en-us/library/ttcz0bys(VS.80).aspx
#  pragma warning(disable:4996)

   // Disables a warning for a this pointer that is passed to a base class in the constructor
   // initializer list.
#  pragma warning(disable:4355)

   // Disables the warning for ignored C++ exception specifications.
#  pragma warning(disable:4290)

#endif




//=================================================================================================
//
//  INTEL WARNINGS
//
//=================================================================================================

#if defined(__INTEL_COMPILER) || defined(__ICL)

   // Disables a 'deprecated' warning for some standard library functions.
#  pragma warning(disable:1786)

#endif

#endif
