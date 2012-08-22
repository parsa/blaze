//=================================================================================================
/*!
//  \file blaze/system/MPI.h
//  \brief System settings for the MPI parallelization
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

#ifndef _BLAZE_SYSTEM_MPI_H_
#define _BLAZE_SYSTEM_MPI_H_


//=================================================================================================
//
//  MPI MODE CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compilation switch for the MPI parallelization.
// \ingroup mpi
//
// This compilation switch enables/disables the MPI parallelization.
//
// Possible settings for the MPI switch:
//  - Deactivated: \b 0
//  - Activated  : \b 1
//
// Note that changing the setting of the MPI parallel mode requires a recompilation of the
// Blaze library. Also note that this switch is automatically set by the configuration script
// of the Blaze library.
*/
#define BLAZE_MPI_PARALLEL_MODE 0
//*************************************************************************************************

#endif
