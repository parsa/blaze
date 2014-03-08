//=================================================================================================
/*!
//  \file blaze/util/Memory.h
//  \brief Header file for memory allocation and deallocation functionality
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
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/Null.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ALLOCATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Aligned array allocation.
// \ingroup util
//
// \param size The number of elements of the given type to allocate.
// \return Pointer to the first element of the aligned array.
//
// The allocate function provides the functionality to allocate memory based on the alignment
// restrictions of the given data type. For instance, in case the given type is a fundamental,
// built-in data type and in case SSE vectorization is possible, the returned memory is guaranteed
// to be at least 16-byte aligned. In case AVX is active, the memory is even guaranteed to be at
// least 32-byte aligned.
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
   void* tmp( NULL );
   const size_t alignment( AlignmentTrait<T>::value );

   if( alignment >= 8UL ) {
#if defined(_MSC_VER)
      tmp = _aligned_malloc( size*sizeof(T), alignment );
      if( tmp != NULL )
#else
      if( !posix_memalign( &tmp, alignment, size*sizeof(T) ) )
#endif
         return reinterpret_cast<T*>( tmp );
      else throw std::bad_alloc();
   }
   else return ::new T[size];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deallocation of memory.
// \ingroup util
//
// \param address The address of the first element of the array to be deallocated.
// \return void
//
// This function deallocates the given memory that was previously allocated via the allocate()
// function.
*/
template< typename T >
void deallocate( T* address )
{
   const size_t alignment( AlignmentTrait<T>::value );

   if( alignment >= 8UL ) {
#if defined(_MSC_VER)
      _aligned_free( address );
#else
      free( address );
#endif
   }
   else delete[] address;
}
//*************************************************************************************************

} // namespace blaze

#endif
