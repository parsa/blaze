//=================================================================================================
/*!
//  \file src/util/Thread.cpp
//  \brief Source file for the Thread class
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

#include <boost/bind.hpp>
#include <blaze/util/Assert.h>
#include <blaze/util/Thread.h>
#include <blaze/util/ThreadPool.h>


namespace blaze {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Starting a thread in a thread pool.
//
// \param pool Handle to the managing thread pool.
//
// This function creates a new thread in the given thread pool. The thread is kept alive until
// explicitly killed by the managing thread pool.
*/
Thread::Thread( ThreadPool* pool )
   : terminated_( false )  // Thread termination flag
   , pool_      (  pool )  // Handle to the managing thread pool
   , thread_    (   0   )  // Handle to the thread of execution
{
   thread_.reset( new boost::thread( boost::bind( &Thread::run, this ) ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Thread class.
*/
Thread::~Thread()
{}
//*************************************************************************************************




//=================================================================================================
//
//  THREAD EXECUTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Execution function for threads in a thread pool.
//
// This function is executed by any thread managed by a thread pool.
*/
void Thread::run()
{
   // Checking the thread pool handle
   BLAZE_INTERNAL_ASSERT( pool_, "Uninitialized pool handle detected" );

   // Executing scheduled tasks
   while( pool_->executeTask() ) {}

   // Setting the termination flag
   terminated_ = true;
}
//*************************************************************************************************

} // namespace blaze
