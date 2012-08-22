//=================================================================================================
/*!
//  \file blaze/config/CacheSize.h
//  \brief Configuration of the cache size of the target architecture
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


namespace blaze {

//*************************************************************************************************
/*!\brief Cache size of the target architecture.
// \ingroup config
//
// This setting specifies the available cache size in Byte of the used target architecture.
// Several algorithms use this setting for an optimized evaluation.
//
// The size of the cache is specified in Byte. For instance, a cache of 3 MiByte must therefore
// be specified as 3145728.
*/
const size_t cacheSize = 3145728UL;
//*************************************************************************************************

} // namespace blaze
