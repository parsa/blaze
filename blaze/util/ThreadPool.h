//=================================================================================================
/*!
//  \file blaze/util/ThreadPool.h
//  \brief Header file of the ThreadPool class
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

#ifndef _BLAZE_UTIL_THREADPOOL_THREADPOOL_H_
#define _BLAZE_UTIL_THREADPOOL_THREADPOOL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/bind.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/mutex.hpp>
#include <blaze/util/NonCopyable.h>
#include <blaze/util/PtrVector.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Thread.h>
#include <blaze/util/threadpool/Task.h>
#include <blaze/util/threadpool/TaskQueue.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup threads Thread parallelization
// \ingroup util
//
// The Blaze library offers the capability to either create individual threads for specific
// tasks (see the Thread class desciption for details) and to create a thread pool according
// to the thread pool pattern (see the ThreadPool class description). Both class descriptions
// offer examples for the setup of threads and the parallel execution of concurrent tasks.
*/
/*!\brief Implementation of a thread pool.
// \ingroup threads
//
// \section threadpool_general General
//
// The ThreadPool class represents a thread pool according to the thread pool pattern (see for
// example http://en.wikipedia.org/wiki/Thread_pool_pattern). It manages a certain number of
// threads in order to process a larger number of independent tasks.
//
// \image html threadpool.png
// \image latex threadpool.eps "Thread pool pattern" width=430pt
//
// The primary purpose of a thread pool is the reuse of system resources: instead of creating
// a single thread for every individual task, threads are reused to handle several tasks. This
// increases the performance in comparison to different threading strategies, as illustrated
// in the graph below. The first bar indicates the sequential performance of 1000 matrix-matrix
// multiplications of arbitrarily sized square matrices. The second bar shows the performance
// of the same work performed by 1000 distinct threads (i.e. one thread for each matrix-matrix
// multiplication) on a quad-core system. In this case, all cores of the system can be used,
// but the additional overhead of creating and managing new threads prevents the expected
// performance increase by a factor of four. The third bar illustrates the performance of four
// threads distributing the work between them (i.e. 250 matrix-matrix multiplications per
// thread), again using the same quad-core system. This approach nearly achieves four times
// the performance of the sequential execution. The fourth bar represents the performance of
// the ThreadPool class using fourth threads for the execution of the 1000 individual
// multiplications.
//
// \image html threadpool2.png
// \image latex threadpool2.eps "Performance comparison of different thread strategies" width=450pt
//
// Additionally, the thread pool approach simplifies load balancing and increases the stability
// of the system.
//
//
// \section threadpool_setup Using the ThreadPool class
//
// The following example demonstrates the use of the ThreadPool class. In contrast to the setup
// of individual threads (see the Thread class description for more details), it is not necessary
// to create and manage individual threads, but only to schedules tasks for the accordingly sized
// thread pool.

   \code
   // Definition of a function with no arguments that returns void
   void function0() { ... }

   // Definition of a functor (function object) taking two arguments and returning void
   struct Functor2
   {
      void operator()( int a, int b ) { ... }
   };

   int main()
   {
      // Creating a thread pool with initially two working threads
      ThreadPool threadpool( 2 );

      // Scheduling two concurrent tasks
      threadpool.schedule( function0 );
      threadpool.schedule( Functor2(), 4, 6 );

      // Waiting for the thread pool to complete both tasks
      threadpool.wait();

      // Resizing the thread pool to four working threads
      threadpool.resize( 4 );

      // Scheduling other concurrent tasks
      ...
      threadpool.schedule( function0 );
      ...

      // At the end of the thread pool scope, all tasks remaining in the task queue are removed
      // and all currently running tasks are completed. Additionally, all acquired resources are
      //safely released.
   }
   \endcode

// Note that the ThreadPool class schedule() function allows for up to five arguments for the
// given functions/functors.
//
//
// \section thread_exception Throwing exceptions in a thread parallel environment
//
// It can happen that during the execution of a given task a thread encounters an erroneous
// situation and has to throw an exception. However, exceptions thrown in the usual way
// cannot be caught by a try-catch-block in the main thread of execution:

   \code
   // Definition of a function throwing a std::runtime_error during its execution
   void task()
   {
      ...
      throw std::runtime_error( ... );
      ...
   }

   // Creating a thread pool executing the throwing function. Although the setup, the scheduling
   // of the task, the wait() function and the destruction of the thread pool are encapsuled
   // inside a try-catch-block, the exception cannot be caught and results in an abortion of the
   // program.
   try {
      Threadpool threadpool( 2 );
      thread.schedule( task );
      threadpool.wait();
   }
   catch( ... )
   {
      ...
   }
   \endcode

// The only possible way to transport exceptions between threads is to use the according boost
// functionality demonstrated in the following example. Note that any function/functor scheduled
// for execution is responsible to handle exceptions in this way!

   \code
   #include <boost/bind.hpp>
   #include <boost/exception_ptr.hpp>

   // Definition of a function that happens to throw an exception. In order to throw the
   // exception, boost::enable_current_exception() is used in combination with throw.
   void throwException()
   {
      ...
      throw boost::enable_current_exception( std::runtime_error( ... ) );
      ...
   }

   // Definition of a thread function. The try-catch-block catches the exception and uses the
   // boost::current_exception() function to get a boost::exception_ptr object.
   void task( boost::exception_ptr& error )
   {
      try {
         throwException();
         error = boost::exception_ptr();
      }
      catch( ... ) {
         error = boost::current_exception();
      }
   }

   // The function that start a thread of execution can pass along a boost::exception_ptr object
   // that is set in case of an exception. Note that boost::current_exception() captures the
   // original type of the exception object. The exception can be thrown again using the
   // boost::rethrow_exception() function.
   void work()
   {
      boost::exception_ptr error;

      ThreadPool threadpool( 2 );
      threadpool.schedule( boost::bind( task, boost::ref(error) ) );
      threadpool.wait();

      if( error ) {
         std::cerr << " Exception during thread execution!\n\n";
         boost::rethrow_exception( error );
      }
   }
   \endcode
*/
class ThreadPool : private NonCopyable
{
 private:
   //**Type definitions****************************************************************************
   typedef PtrVector<Thread>          Threads;    //!< Type of the thread container.
   typedef threadpool::TaskQueue      TaskQueue;  //!< Type of the task queue.
   typedef boost::mutex               Mutex;      //!< Type of the mutex.
   typedef Mutex::scoped_lock         Lock;       //!< Type of a locking object.
   typedef boost::condition_variable  Condition;  //!< Condition variable type.
   //**********************************************************************************************

 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit ThreadPool( size_t n );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ThreadPool();
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline bool   isEmpty()  const;
   inline size_t size()     const;
   inline size_t active()   const;
   inline size_t ready()    const;
   //@}
   //**********************************************************************************************

   //**Scheduling functions************************************************************************
   /*!\name Scheduling functions */
   //@{
   template< typename Callable >
   void schedule( Callable func );

   template< typename Callable, typename A1 >
   void schedule( Callable func, A1 a1 );

   template< typename Callable, typename A1, typename A2 >
   void schedule( Callable func, A1 a1, A2 a2 );

   template< typename Callable, typename A1, typename A2, typename A3 >
   void schedule( Callable func, A1 a1, A2 a2, A3 a3 );

   template< typename Callable, typename A1, typename A2, typename A3, typename A4 >
   void schedule( Callable func, A1 a1, A2 a2, A3 a3, A4 a4 );

   template< typename Callable, typename A1, typename A2, typename A3, typename A4, typename A5 >
   void schedule( Callable func, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5 );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void resize( size_t n );
   void wait();
   void clear();
   //@}
   //**********************************************************************************************

 private:
   //**Thread functions****************************************************************************
   /*!\name Thread functions */
   //@{
   void createThread();
   bool executeTask();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   volatile size_t total_;     //!< Total number of threads in the thread pool.
   volatile size_t expected_;  //!< Expected number of threads in the thread pool.
                               /*!< This number may differ from the total number of threads
                                    during a resize of the thread pool. */
   volatile size_t active_;    //!< Number of currently active/busy threads.
   Threads threads_;           //!< The threads contained in the thread pool.
   TaskQueue taskqueue_;       //!< Task queue for the scheduled tasks.
   mutable Mutex mutex_;       //!< Synchronization mutex.
   Condition waitForTask_;     //!< Wait condition for idle threads.
   Condition waitForThread_;   //!< Wait condition for the thread management.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   friend class Thread;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether any tasks are scheduled for execution.
//
// \return \a true in case task are scheduled, \a false otherwise.
*/
inline bool ThreadPool::isEmpty() const
{
   Lock lock( mutex_ );
   return taskqueue_.isEmpty();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size of the thread pool.
//
// \return The total number of threads in the thread pool.
*/
inline size_t ThreadPool::size() const
{
   Lock lock( mutex_ );
   return expected_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of currently active/busy threads.
//
// \return The number of currently active threads.
*/
inline size_t ThreadPool::active() const
{
   Lock lock( mutex_ );
   return active_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of currently ready/inactive threads.
//
// \return The number of currently ready threads.
*/
inline size_t ThreadPool::ready() const
{
   Lock lock( mutex_ );
   return expected_ - active_;
}
//*************************************************************************************************




//=================================================================================================
//
//  SCHEDULING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scheduling the given zero argument function/functor for execution.
//
// \param func The given function/functor.
// \return void
//
// This function schedules the given function/functor for execution. The given function/functor
// must be copyable, must be callable without arguments and must return void.
*/
template< typename Callable >  // Task type
void ThreadPool::schedule( Callable func )
{
   Lock lock( mutex_ );
   taskqueue_.push( func );
   waitForTask_.notify_one();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scheduling the given unary function/functor for execution.
//
// \param func The given function/functor.
// \param a1 The first argument.
// \return void
//
// This function schedules the given function/functor for execution. The given function/functor
// must be copyable, must be callable with one argument and must return void.
*/
template< typename Callable  // Type of the function/functor
        , typename A1 >      // Type of the first argument
void ThreadPool::schedule( Callable func, A1 a1 )
{
   Lock lock( mutex_ );
   taskqueue_.push( boost::bind<void>( func, a1 ) );
   waitForTask_.notify_one();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scheduling the given binary function/functor for execution.
//
// \param func The given function/functor.
// \param a1 The first argument.
// \param a2 The second argument.
// \return void
//
// This function schedules the given function/functor for execution. The given function/functor
// must be copyable, must be callable with two arguments and must return void.
*/
template< typename Callable  // Type of the function/functor
        , typename A1        // Type of the first argument
        , typename A2 >      // Type of the second argument
void ThreadPool::schedule( Callable func, A1 a1, A2 a2 )
{
   Lock lock( mutex_ );
   taskqueue_.push( boost::bind<void>( func, a1, a2 ) );
   waitForTask_.notify_one();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scheduling the given ternary function/functor for execution.
//
// \param func The given function/functor.
// \param a1 The first argument.
// \param a2 The second argument.
// \param a3 The third argument.
// \return void
//
// This function schedules the given function/functor for execution. The given function/functor
// must be copyable, must be callable with three arguments and must return void.
*/
template< typename Callable  // Type of the function/functor
        , typename A1        // Type of the first argument
        , typename A2        // Type of the second argument
        , typename A3 >      // Type of the third argument
void ThreadPool::schedule( Callable func, A1 a1, A2 a2, A3 a3 )
{
   Lock lock( mutex_ );
   taskqueue_.push( boost::bind<void>( func, a1, a2, a3 ) );
   waitForTask_.notify_one();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scheduling the given four argument function/functor for execution.
//
// \param func The given function/functor.
// \param a1 The first argument.
// \param a2 The second argument.
// \param a3 The third argument.
// \param a4 The fourth argument.
// \return void
//
// This function schedules the given function/functor for execution. The given function/functor
// must be copyable, must be callable with four arguments and must return void.
*/
template< typename Callable  // Type of the function/functor
        , typename A1        // Type of the first argument
        , typename A2        // Type of the second argument
        , typename A3        // Type of the third argument
        , typename A4 >      // Type of the fourth argument
void ThreadPool::schedule( Callable func, A1 a1, A2 a2, A3 a3, A4 a4 )
{
   Lock lock( mutex_ );
   taskqueue_.push( boost::bind<void>( func, a1, a2, a3, a4 ) );
   waitForTask_.notify_one();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scheduling the given four argument function/functor for execution.
//
// \param func The given function/functor.
// \param a1 The first argument.
// \param a2 The second argument.
// \param a3 The third argument.
// \param a4 The fourth argument.
// \param a5 The fifth argument.
// \return void
//
// This function schedules the given function/functor for execution. The given function/functor
// must be copyable, must be callable with five arguments and must return void.
*/
template< typename Callable  // Type of the function/functor
        , typename A1        // Type of the first argument
        , typename A2        // Type of the second argument
        , typename A3        // Type of the third argument
        , typename A4        // Type of the fourth argument
        , typename A5 >      // Type of the fifth argument
void ThreadPool::schedule( Callable func, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5 )
{
   Lock lock( mutex_ );
   taskqueue_.push( boost::bind<void>( func, a1, a2, a3, a4, a5 ) );
   waitForTask_.notify_one();
}
//*************************************************************************************************

} // namespace blaze

#endif
