//=================================================================================================
/*!
//  \file src/util/SystemClock.cpp
//  \brief Source file for the SystemClock class
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
// Platform/compiler-specific includes
//*************************************************************************************************

#include <blaze/system/WarningDisable.h>


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Null.h>
#include <blaze/util/SystemClock.h>


namespace blaze {

//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

time_t SystemClock::start_( std::time( NULL ) );




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SystemClock class.
*/
SystemClock::SystemClock()
   : Singleton<SystemClock>()  // Initialization of the Singleton base object
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the SystemClock class.
*/
SystemClock::~SystemClock()
{}
//*************************************************************************************************

} // namespace blaze
