//=================================================================================================
/*!
//  \file blaze/util/AlignedStorage.h
//  \brief Header file for the AlignedStorage implementation
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

#ifndef _BLAZE_UTIL_ALIGNEDSTORAGE_H_
#define _BLAZE_UTIL_ALIGNEDSTORAGE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Alignment.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/typetraits/AlignmentOf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the AlignedStorage class template.
// \ingroup util
*/
template< size_t Alignment >
struct AlignedStorageHelper;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 1-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 1UL ) AlignedStorageHelper<1UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 2-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 2UL ) AlignedStorageHelper<2UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 4-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 4UL ) AlignedStorageHelper<4UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 8-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 8UL ) AlignedStorageHelper<8UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 16-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 16UL ) AlignedStorageHelper<16UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 32-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 32UL ) AlignedStorageHelper<32UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 64-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 64UL ) AlignedStorageHelper<64UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 128-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 128UL ) AlignedStorageHelper<128UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 256-byte alignment.
// \ingroup util
*/
template<>
struct BLAZE_ALIGN( 256UL ) AlignedStorageHelper<256UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief POD data type with a fixed alignment.
// \ingroup util
//
// The AlignedStorage class template represents a POD data type with a fixed alignment. Via this
// class it is possible to enforce a specific, type-based alignment for static data types. The
// required alignment is evaluated based on the given data type \a T. In case \a T is a built-in,
// vectorizable data type, AlignedStorage enforces an alignment of 16 or 32 bytes, depending on
// the active SSE/AVX level. In all other cases, no specific alignment is enforced.
//
// In the following code example, the StaticVector class, representing a vector of N statically
// allocated elements, is non-publicly deriving from the AlignedStorage class template. By this,
// the N data elements of type \a T are aligned according to the requirements of \a T.

   \code
   template< typename T, size_t N >
   class StaticVector : private AlignedStorage<T>
   {
      // ...
      T v_[N];
   };
   \endcode

// Note that since the AlignedStorage class template is designed as a base class proper alignment
// is heavily depending on the compiler's memory model for base classes (which influences how the
// base class's alignment restrictions are transfered to the deriving class) and it's capacity to
// perform an empty base class optimization (EBO). Currently, this approach works with the GCC,
// the Intel, and the Visual Studio compiler, but fails with the Clang compiler.
*/
template< typename T >
class AlignedStorage : private AlignedStorageHelper< AlignmentOf<T>::value >
{
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
