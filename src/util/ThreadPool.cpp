//=================================================================================================
/*!
//  \file src/util/ThreadPool.cpp
//  \brief Source file of the ThreadPool class
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

#include <stdexcept>
#include <boost/scoped_ptr.hpp>
#include <blaze/util/Assert.h>
#include <blaze/util/ThreadPool.h>


namespace blaze {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the ThreadPool class.
//
// \param n Initial number of threads \f$[1..\infty)\f$.
//
// This constructor creates a thread pool with initially \a n new threads. All threads are
// initially idle until a task is scheduled.
*/
ThreadPool::ThreadPool( size_t n )
   : total_     ( 0 )  // Total number of threads in the thread pool
   , expected_  ( 0 )  // Expected number of threads in the thread pool
   , active_    ( 0 )  // Number of currently active/busy threads
   , threads_      ()  // The threads contained in the thread pool
   , taskqueue_    ()  // Task queue for the scheduled tasks
   , mutex_        ()  // Synchronization mutex
   , waitForTask_  ()  // Wait condition for idle threads
   , waitForThread_()  // Wait condition for the thread management
{
   resize( n );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the ThreadPool class.
//
// The destructor clears all remaining tasks from the task queue and waits for the currently
// active threads to complete their tasks.
*/
ThreadPool::~ThreadPool()
{
   Lock lock( mutex_ );

   // Removing all currently queued tasks
   taskqueue_.clear();

   // Setting the expected number of threads
   expected_ = 0;

   // Notifying all idle threads
   waitForTask_.notify_all();

   // Waiting for all threads to terminate
   while( total_ != 0 ) {
      waitForThread_.wait( lock );
   }

   // Joining all threads
   typedef Threads::Iterator  Iterator;
   const Iterator end( threads_.end() );
   for( Iterator thread=threads_.begin(); thread!=end; ++thread ) {
      thread->join();
   }

   // Destroying all threads
   threads_.clear();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Changes the total number of threads in the thread pool.
//
// \param n The new number of threads \f$[1..\infty)\f$.
// \return void
// \exception std::invalid_argument Invalid number of threads.
//
// This function changes the size of the thread pool, i.e. changes the total number of threads
// contained in the pool. If \a n is smaller than the current size of the thread pool, the
// according number of threads is removed from the pool, otherwise new threads are added to
// the pool.
*/
void ThreadPool::resize( size_t n )
{
   // Checking the given number of threads
   if( n == 0 )
      throw std::invalid_argument( "Invalid number of threads" );

   // Adjusting the number of threads
   {
      Lock lock( mutex_ );

      // Adding new threads to the thread pool
      if( n > expected_ ) {
         for( size_t i=expected_; i<n; ++i )
            createThread();
      }

      // Removing threads from the pool
      else {
         expected_ = n;
         waitForTask_.notify_all();
      }
   }

   // Joining and destroying any terminated thread
   typedef Threads::Iterator  Iterator;
   for( Iterator thread=threads_.begin(); thread!=threads_.end(); ) {
      if( thread->hasTerminated() ) {
         thread->join();
         thread = threads_.erase( thread );
      }
      else ++thread;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Waiting for all scheduled tasks to be completed.
//
// \return void
//
// This function blocks until all scheduled tasks have been completed.
*/
void ThreadPool::wait()
{
   Lock lock( mutex_ );

   while( !taskqueue_.isEmpty() || active_ > 0 ) {
      waitForThread_.wait( lock );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all scheduled tasks from the thread pool.
//
// \return void
//
// This function removes all currently scheduled tasks from the thread pool. The total number
// of threads remains unchanged and all active threads continue completing their tasks.
*/
void ThreadPool::clear()
{
   Lock lock( mutex_ );
   taskqueue_.clear();
}
//*************************************************************************************************




//=================================================================================================
//
//  THREAD FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adding a new thread to the thread pool.
//
// \return void
*/
void ThreadPool::createThread()
{
   threads_.pushBack( new Thread( this ) );
   ++total_;
   ++expected_;
   ++active_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Executing a scheduled task.
//
// \return void
//
// This function is repeatedly called by every thread to execute one of the scheduled tasks.
// In case there is no task available, the thread blocks and waits for a new task to be
// scheduled.
*/
bool ThreadPool::executeTask()
{
   threadpool::Task task;

   // Acquiring a scheduled task
   {
      Lock lock( mutex_ );

      while( taskqueue_.isEmpty() )
      {
         --active_;
         waitForThread_.notify_all();

         if( total_ > expected_ ) {
            --total_;
            return false;
         }

         waitForTask_.wait( lock );
         ++active_;
      }

      BLAZE_INTERNAL_ASSERT( !taskqueue_.isEmpty(), "Empty task queue detected" );
      task = taskqueue_.pop();
   }

   // Executing the task
   task();

   return true;
}
//*************************************************************************************************

} // namespace blaze
