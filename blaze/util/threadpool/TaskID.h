//=================================================================================================
/*!
//  \file blaze/util/threadpool/TaskID.h
//  \brief Implementation of a task handle
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

#ifndef _BLAZE_UTIL_THREADPOOL_TASKID_H_
#define _BLAZE_UTIL_THREADPOOL_TASKID_H_


namespace blaze {

namespace threadpool {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

class Task;




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Handle for a task object.
// \ingroup threads
*/
typedef Task*  TaskID;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Handle for a constant task object.
// \ingroup threads
*/
typedef const Task*  ConstTaskID;
//*************************************************************************************************

} // namespace threadpool

} // namespace blaze

#endif
