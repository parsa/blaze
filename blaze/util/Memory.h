//=================================================================================================
/*!
//  \file blaze/util/Memory.h
//  \brief Header file for memory allocation and deallocation functionality
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

#ifndef _BLAZE_UTIL_MEMORY_H_
#define _BLAZE_UTIL_MEMORY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#if defined(_MSC_VER)
#  include <malloc.h>
#endif
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <boost/checked_delete.hpp>
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/Null.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsVectorizable.h>


namespace blaze {

//=================================================================================================
//
//  ALLOCATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Aligned array allocation.
//
// \param size The number of elements of the given type to allocate.
// \return Pointer to the first element of the aligned array.
//
// The allocate function provides the functionality to allocate memory based on the alignment
// restrictions of the given data type. In case the given type is a fundamental, built-in data
// type and in case SSE vectorization is possible, the returned memory is guaranteed to be at
// least 16-byte aligned. In case AVX in active, the memory is even guaranteed to be 32-byte
// aligned. For all other, non-builtin data types, the system-specific alignment strategy is
// used.
//
// Examples:

   \code
   // Guaranteed to be 16-byte aligned (32-byte aligned in case AVX is used)
   double* dp = allocate<double>( 10UL );
   \endcode
*/
template< typename T >
T* allocate( size_t size )
{
   if( IsVectorizable<T>::value )
   {
      void* tmp( NULL );
      const size_t alignment( AlignmentTrait<T>::value );

#if defined(_MSC_VER)
      tmp = _aligned_malloc( size*sizeof(T), alignment );
      if( tmp != NULL )
#else
      if( !posix_memalign( &tmp, alignment, size*sizeof(T) ) )
#endif
         return reinterpret_cast<T*>( tmp );
      else throw std::bad_alloc();
   }
   else return new T[size];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deallocation of memory.
//
// \param address The address of the first element of the array to be deallocated.
// \return void
*/
template< typename T >
void deallocate( T* address )
{
   if( IsVectorizable<T>::value && address != NULL ) {
#if defined(_MSC_VER)
      _aligned_free( address );
#else
      free( address );
#endif
   }
   else {
      boost::checked_array_delete( address );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
