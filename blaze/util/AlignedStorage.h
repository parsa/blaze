//=================================================================================================
/*!
//  \file blaze/util/AlignedStorage.h
//  \brief Header file for the AlignedStorage implementation
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

#ifndef _BLAZE_UTIL_ALIGNEDSTORAGE_H_
#define _BLAZE_UTIL_ALIGNEDSTORAGE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Alignment.h>
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Volatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the AlignedStorage class template.
// \ingroup intrinsics
*/
template< size_t Alignment >
struct AlignedStorageHelper;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 1-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 1UL ) AlignedStorageHelper<1UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 2-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 2UL ) AlignedStorageHelper<2UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 4-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 4UL ) AlignedStorageHelper<4UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 8-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 8UL ) AlignedStorageHelper<8UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 16-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 16UL ) AlignedStorageHelper<16UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 32-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 32UL ) AlignedStorageHelper<32UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedStorage for 64-byte alignment.
// \ingroup intrinsics
*/
template<>
struct BLAZE_ALIGN( 64UL ) AlignedStorageHelper<64UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief POD data type with a fixed alignment.
// \ingroup intrinsics
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
*/
template< typename T >
class AlignedStorage : private AlignedStorageHelper< AlignmentTrait<T>::value >
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
