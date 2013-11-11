//=================================================================================================
/*!
//  \file src/util/Thread.cpp
//  \brief Source file for the Thread class
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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
