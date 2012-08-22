//=================================================================================================
/*!
//  \file blazemark/system/Armadillo.h
//  \brief System settings for the Armadillo performance benchmarks
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

#ifndef _BLAZEMARK_SYSTEM_ARMADILLO_H_
#define _BLAZEMARK_SYSTEM_ARMADILLO_H_


//=================================================================================================
//
//  ARMADILLO CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compilation switch for the Armadillo performance benchmarks.
//
// This compilation switch enables/disables the Armadillo library performance benchmarks.
//
// Possible settings for the Armadillo switch:
//  - Deactivated: \b 0
//  - Activated  : \b 1
//
// Note that changing this setting requires a recompilation of the Blaze benchmark suite. Also
// note that this switch is automatically set by the configuration script of the Blaze benchmark
// suite.
*/
#define BLAZEMARK_ARMADILLO_MODE 0
//*************************************************************************************************

#endif
