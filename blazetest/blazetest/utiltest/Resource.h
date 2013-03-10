//=================================================================================================
/*!
//  \file blazetest/utiltest/Resource.h
//  \brief Header file for Resource class template
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

#ifndef _BLAZETEST_UTILTEST_RESOURCE_H_
#define _BLAZETEST_UTILTEST_RESOURCE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blazetest/utiltest/InstanceCounter.h>


namespace blazetest {

namespace utiltest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of an instance counted resource.
//
// The Resource class represents an important resource for testing purposes. It is instance
// counted via the InstanceCounter class.
*/
class Resource : public InstanceCounter<Resource>
{};
//*************************************************************************************************

} // namespace utiltest

} // namespace blazetest

#endif
