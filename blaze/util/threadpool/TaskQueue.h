//=================================================================================================
/*!
//  \file blaze/util/threadpool/TaskQueue.h
//  \brief Task queue for the thread pool
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

#ifndef _BLAZE_UTIL_THREADPOOL_TASKQUEUE_H_
#define _BLAZE_UTIL_THREADPOOL_TASKQUEUE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <deque>
#include <blaze/util/policies/PtrDelete.h>
#include <blaze/util/threadpool/Task.h>
#include <blaze/util/threadpool/TaskID.h>
#include <blaze/util/UniquePtr.h>


namespace blaze {

namespace threadpool {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Task queue for the thread pool.
// \ingroup threads
//
// The TaskQueue class represents the internal task container of a thread pool. It uses a FIFO
// (first in, first out) strategy to store and remove the assigned tasks.
*/
class TaskQueue
{
 private:
   //**Type definitions****************************************************************************
   typedef std::deque<Task*>  Tasks;  //!< FIFO container for tasks.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef Tasks::size_type  SizeType;  //!< Size type of the task queue.
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline TaskQueue();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~TaskQueue();
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline SizeType maxSize()  const;
   inline SizeType size()     const;
   inline bool     isEmpty()  const;
   //@}
   //**********************************************************************************************

   //**Element functions***************************************************************************
   /*!\name Element functions */
   //@{
   inline void   push ( TaskID task );
   inline TaskID pop  ();
   inline void   clear();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void swap( TaskQueue& tq ) /* throw() */;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Tasks tasks_;  //!< FIFO container for the contained tasks.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for TaskQueue.
*/
inline TaskQueue::TaskQueue()
   : tasks_()  // FIFO container for the contained tasks
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the TaskQueue class.
//
// The destructor destroys any remaining task in the task queue.
*/
TaskQueue::~TaskQueue()
{
   clear();
}
//*************************************************************************************************




//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the maximum possible size of a task queue.
//
// \return The maximum possible size.
*/
inline TaskQueue::SizeType TaskQueue::maxSize() const
{
   return tasks_.max_size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size of the task queue.
//
// \return The current size.
//
// This function returns the number of the currently contained tasks.
*/
inline TaskQueue::SizeType TaskQueue::size() const
{
   return tasks_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns \a true if the task queue has no elements.
//
// \return \a true if the task queue is empty, \a false if it is not.
*/
inline bool TaskQueue::isEmpty() const
{
   return tasks_.empty();
}
//*************************************************************************************************




//=================================================================================================
//
//  ELEMENT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adding a task to the end of the task queue.
//
// \param task The task to be added to the end of the task queue.
// \return void
//
// This function adds the given task to the end of the task queue. It runs in constant time.
*/
inline void TaskQueue::push( TaskID task )
{
   UniquePtr<Task> ptr( task );
   tasks_.push_back( task );
   ptr.release();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the task from the front of the task queue.
//
// \return The first task in the task queue.
*/
inline TaskID TaskQueue::pop()
{
   TaskID task( tasks_.front() );
   tasks_.pop_front();
   return task;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all tasks from the task queue.
//
// \return void
*/
inline void TaskQueue::clear()
{
   std::for_each( tasks_.begin(), tasks_.end(), PtrDelete() );
   tasks_.clear();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Swapping the contents of two task queues.
//
// \param tq The task queue to be swapped.
// \return void
// \exception no-throw guarantee.
*/
inline void TaskQueue::swap( TaskQueue& tq ) /* throw() */
{
   tasks_.swap( tq.tasks_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name TaskQueue operators */
//@{
inline void swap( TaskQueue& a, TaskQueue& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two task queues.
//
// \param a The first task queue to be swapped.
// \param b The second task queue to be swapped.
// \return void
// \exception no-throw guarantee.
*/
inline void swap( TaskQueue& a, TaskQueue& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************

} // namespace threadpool

} // namespace blaze


#endif
