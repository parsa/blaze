//=================================================================================================
/*!
//  \file blaze/util/AlignedArray.h
//  \brief Header file for the AlignedArray implementation
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

#ifndef _BLAZE_UTIL_ALIGNEDARRAY_H_
#define _BLAZE_UTIL_ALIGNEDARRAY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Alignment.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
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
/*!\brief Auxiliary helper struct for the AlignedArray class template.
// \ingroup util
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
struct AlignedArrayHelper;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 1-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,1UL>
{
   BLAZE_ALIGN( 1UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 2-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,2UL>
{
   BLAZE_ALIGN( 2UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 4-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,4UL>
{
   BLAZE_ALIGN( 4UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 8-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,8UL>
{
   BLAZE_ALIGN( 8UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 16-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,16UL>
{
   BLAZE_ALIGN( 16UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 32-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,32UL>
{
   BLAZE_ALIGN( 32UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 64-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,64UL>
{
   BLAZE_ALIGN( 64UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 128-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,128UL>
{
   BLAZE_ALIGN( 128UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignedArrayHelper for 256-byte alignment.
// \ingroup util
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
struct AlignedArrayHelper<Type,N,256UL>
{
   BLAZE_ALIGN( 256UL ) Type v_[N];  //!< The statically allocated array elements.
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of a static array with a fixed alignment.
// \ingroup util
//
// The AlignedArray class template represents a static array with a guaranteed, fixed alignment.
// The type of the array elements, the number of elements and the alignment of the array can be
// specified via the three template parameters:

   \code
   template< typename Type, size_t N, size_t Alignment >
   class AlignedArray;
   \endcode

// The alignment of the array, which must be a power of two (i.e. 1, 2, 4, 8, ...), can either be
// specified explicitly via the template parameter \a Alignment or it is evaluated automatically
// based on the alignment requirements of the given data type \a Type. In the latter case, if
// \a T is a built-in, vectorizable data type, AlignedArray enforces an alignment of 16 or 32
// bytes, depending on the active SSE/AVX level. In all other cases, no specific alignment is
// enforced.
//
// AlignedArray can be used exactly like any built-in static array. It is possible to access the
// individual element via the subscript operator and the array can be used wherever a pointer is
// expected:

   \code
   void func( const int* );

   blaze::AlignedArray<int,100UL> array;

   array[10] = 2;  // Accessing and assigning the 10th array element
   func( array );  // Passing the aligned array to a function expecting a pointer
   \endcode

// Currently, the only limitation is that an aligned array cannot be statically initialized via
// array initialization:

   \code
   blaze::AlignedArray<int,3UL> array = { 1, 2, 3 };  // Compilation error!
   \endcode

// This limitation will be resolved when AlignedArray is updated to C++11.
*/
template< typename Type                                  // Data type of the elements
        , size_t N                                       // Number of elements
        , size_t Alignment = AlignmentOf<Type>::value >  // Array alignment
class AlignedArray
{
 public:
   //**Type definitions****************************************************************************
   typedef Type*        Pointer;         //!< Pointer to a non-constant array element.
   typedef const Type*  ConstPointer;    //!< Pointer to a constant array element.
   typedef Type&        Reference;       //!< Reference to a non-constant array element.
   typedef const Type&  ConstReference;  //!< Reference to a constant array element.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline AlignedArray();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator Pointer     ();
   inline operator ConstPointer() const;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Pointer        data();
   inline ConstPointer   data() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\name Member variables */
   //@{
   AlignedArrayHelper<Type,N,Alignment> array_;  //!< The aligned array of size N.
   //@}
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for AlignedArray.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline AlignedArray<Type,N,Alignment>::AlignedArray()
{
   BLAZE_INTERNAL_ASSERT( checkAlignment( array_.v_ ), "Invalid alignment detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to a pointer.
//
// \return The raw pointer of the aligned array.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline AlignedArray<Type,N,Alignment>::operator Pointer()
{
   return array_.v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to a pointer-to-const.
//
// \return The raw pointer of the aligned array.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline AlignedArray<Type,N,Alignment>::operator ConstPointer() const
{
   return array_.v_;
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the array elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// In case BLAZE_USER_ASSERT() is active, this operator performs an index check.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline typename AlignedArray<Type,N,Alignment>::Reference
   AlignedArray<Type,N,Alignment>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < N, "Invalid array access index" );
   return array_.v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the array elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference-to-const to the accessed value.
//
// In case BLAZE_USER_ASSERT() is active, this operator performs an index check.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline typename AlignedArray<Type,N,Alignment>::ConstReference
   AlignedArray<Type,N,Alignment>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < N, "Invalid array access index" );
   return array_.v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the array elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the aligned array.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline typename AlignedArray<Type,N,Alignment>::Pointer
   AlignedArray<Type,N,Alignment>::data()
{
   return array_.v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the array elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the aligned array.
*/
template< typename Type       // Data type of the elements
        , size_t N            // Number of elements
        , size_t Alignment >  // Array alignment
inline typename AlignedArray<Type,N,Alignment>::ConstPointer
   AlignedArray<Type,N,Alignment>::data() const
{
   return array_.v_;
}
//*************************************************************************************************

} // namespace blaze

#endif
