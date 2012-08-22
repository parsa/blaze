//=================================================================================================
/*!
//  \file blaze/util/threadpool/Task.h
//  \brief Header file for the Task base class
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

#ifndef _BLAZE_UTIL_THREADPOOL_TASK_H_
#define _BLAZE_UTIL_THREADPOOL_TASK_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/NonCopyable.h>


namespace blaze {

namespace threadpool {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for executable user tasks.
// \ingroup threads
//
// The Task class represents the base class for all user tasks.
*/
class Task : private NonCopyable
{
 public:
   //**Constructor*********************************************************************************
   // No explicitly declared constructor.
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~Task();
   //@}
   //**********************************************************************************************

   //**Execution functions*************************************************************************
   /*!\name Execution functions */
   //@{
   virtual void run() = 0;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace threadpool

} // namespace blaze

#endif
