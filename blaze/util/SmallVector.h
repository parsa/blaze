//=================================================================================================
/*!
//  \file blaze/util/SmallVector.h
//  \brief Header file for the SmallVector implementation
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_UTIL_SMALLVECTOR_H_
#define _BLAZE_UTIL_SMALLVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <memory>
#include <blaze/util/algorithms/Destroy.h>
#include <blaze/util/algorithms/DestroyAt.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/UninitializedDefaultConstruct.h>
#include <blaze/util/algorithms/UninitializedMove.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/Exception.h>
#include <blaze/util/InitializerList.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blaze/util/typetraits/IsConstructible.h>
#include <blaze/util/typetraits/IsAssignable.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of a dynamic vector with small vector optimization.
// \ingroup util
//
// The SmallVector class template is a hybrid data structure between a static array and a dynamic
// vector. It provides static, in-place memory for up to \a N elements of type \a T, but can grow
// beyond this size by allocating dynamic memory via its allocator of type \a A.
*/
template< typename T                        // Data type of the elements
        , size_t N                          // Number of preallocated elements
        , typename A = std::allocator<T> >  // Type of the allocator
class SmallVector
   : private A
{
 public:
   //**Type definitions****************************************************************************
   using ElementType    = T;         //!< Type of the vector elements.
   using Pointer        = T*;        //!< Pointer to a non-constant vector element.
   using ConstPointer   = const T*;  //!< Pointer to a constant vector element.
   using Reference      = T&;        //!< Reference to a non-constant vector element.
   using ConstReference = const T&;  //!< Reference to a constant vector element.
   using Iterator       = T*;        //!< Iterator over non-constant elements.
   using ConstIterator  = const T*;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SmallVector( const A& alloc = A() );
   explicit inline SmallVector( size_t n, const A& alloc = A() );
   explicit inline SmallVector( size_t n, const T& init, const A& alloc = A() );

   template< typename InputIt >
   explicit inline SmallVector( InputIt first, InputIt last, const A& alloc = A() );

   template< typename U >
   explicit inline SmallVector( initializer_list<U> list, const A& alloc = A() );

   inline SmallVector( const SmallVector& sv );
   inline SmallVector( SmallVector&& sv );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~SmallVector();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index ) noexcept;
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data() noexcept;
   inline ConstPointer   data() const noexcept;
   inline Iterator       begin () noexcept;
   inline ConstIterator  begin () const noexcept;
   inline ConstIterator  cbegin() const noexcept;
   inline Iterator       end   () noexcept;
   inline ConstIterator  end   () const noexcept;
   inline ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   template< typename U >
   inline SmallVector& operator=( initializer_list<U> list );

   inline SmallVector& operator=( const SmallVector& rhs );
   inline SmallVector& operator=( SmallVector&& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool   empty() const noexcept;
   inline size_t size() const noexcept;
   inline size_t capacity() const noexcept;

   inline void     clear();
          void     resize( size_t n );
          void     resize( size_t n, const T& value );
          void     reserve( size_t n );
          void     shrinkToFit();
          void     pushBack( const T& value );
          void     pushBack( T&& value );
          Iterator insert( Iterator pos, const T& value );
          Iterator insert( Iterator pos, T&& value );
          Iterator erase( Iterator pos );
          Iterator erase( Iterator first, Iterator last );
          void     swap( SmallVector& sv ) noexcept( IsNothrowMoveConstructible_v<T> );
   //@}
   //**********************************************************************************************

 private:
   //**Uninitialized struct definition*************************************************************
   /*!\brief Definition of the nested auxiliary struct Uninitialized.
   */
   struct Uninitialized {};
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SmallVector( size_t n, const A& alloc, Uninitialized );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using A::allocate;
   using A::deallocate;

   inline bool isDynamic() const noexcept;
   //@}
   //**********************************************************************************************

   //**********************************************************************************************
   //! Adjustment of the size of the static storage.
   static constexpr size_t NN = ( N > 0UL ? N*sizeof(T) : 1UL );
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   alignas( AlignmentOf_v<T> ) byte_t v_[NN];  //!< The static storage.

   T* begin_;  //!< Pointer to the beginning of the currently used storage.
   T* end_;    //!< Pointer to the end of the currently used storage.
   T* final_;  //!< Pointer to the very end of the currently used storage.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
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
/*!\brief The (default) constructor for SmallVector.
//
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::SmallVector( const A& alloc )
   : SmallVector( 0UL, alloc, Uninitialized() )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a vector of size \a n. No element initialization is performed!
//
// \param n The initial size of the vector.
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::SmallVector( size_t n, const A& alloc )
   : SmallVector( n, alloc, Uninitialized() )
{
   std::uninitialized_fill( begin_, end_, T() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a vector of size \a n.
//
// \param n The initial size of the vector.
// \param init The initial value of the vector elements.
// \param alloc Allocator for all memory allocations of this container.
//
// All vector elements are initialized with the specified value.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::SmallVector( size_t n, const T& init, const A& alloc )
   : SmallVector( n, alloc, Uninitialized() )
{
   std::uninitialized_fill( begin_, end_, init );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a range of elements.
//
// \param first Iterator to the be first element of the input range.
// \param last Iterator to the element one past the last element of the input range.
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T          // Data type of the elements
        , size_t N            // Number of preallocated elements
        , typename A >        // Type of the allocator
template< typename InputIt >  // Type of the iterators
inline SmallVector<T,N,A>::SmallVector( InputIt first, InputIt last, const A& alloc )
   : SmallVector( std::distance( first, last ), alloc, Uninitialized() )
{
   std::uninitialized_copy( first, last, begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all vector elements.
//
// \param list The initializer list.
// \param alloc Allocator for all memory allocations of this container.
//
// This constructor provides the option to explicitly initialize the elements of the small vector
// within a constructor call:

   \code
   blaze::SmallVector<double,8UL> v1{ 4.2, 6.3, -1.2 };
   \endcode

// The vector is sized according to the size of the initializer list and all its elements are
// copy initialized by the elements of the given initializer list.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
template< typename U >  // Type of the initializer list elements
inline SmallVector<T,N,A>::SmallVector( initializer_list<U> list, const A& alloc )
   : SmallVector( list.size(), alloc, Uninitialized() )
{
   std::uninitialized_copy( list.begin(), list.end(), begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for SmallVector.
//
// \param sv The small vector to be copied.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::SmallVector( const SmallVector& sv )
   : SmallVector( sv.size(), A(), Uninitialized() )
{
   std::uninitialized_copy( sv.begin(), sv.end(), begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for SmallVector.
//
// \param sv The small vector to be moved into this instance.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::SmallVector( SmallVector&& sv )
   // arr_ is intentionally not initialized
   : A     ()             // Base class initialization
   , begin_( sv.begin_ )  // Pointer to the beginning of the currently used storage
   , end_  ( sv.end_   )  // Pointer to the end of the currently used storage
   , final_( sv.final_ )  // Pointer to the very end of the currently used storage
{
   if( !sv.isDynamic() ) {
      begin_ = reinterpret_cast<T*>( v_ );
      end_   = begin_ + sv.size();
      final_ = begin_ + N;
      uninitialized_move( sv.begin(), sv.end(), begin_ );
   }

   sv.begin_ = reinterpret_cast<T*>( sv.v_ );
   sv.end_   = sv.begin_;
   sv.final_ = sv.begin_ + N;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary constructor for SmallVector.
//
// \param n The initial size of the vector.
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::SmallVector( size_t n, const A& alloc, Uninitialized )
   // v_ is intentionally not initialized
   // begin_ is intentionally not initialized
   // end_ is intentionally not initialized
   // final_ is intentionally not initialized
   : A( alloc )  // Base class initialization
{
   if( n <= N ) {
      begin_ = reinterpret_cast<T*>( v_ );
      end_   = begin_ + n;
      final_ = begin_ + N;
   }
   else {
      begin_ = allocate( n );
      end_   = begin_ + n;
      final_ = begin_ + n;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for DynamicVector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>::~SmallVector()
{
   using blaze::destroy;

   destroy( begin_, end_ );

   if( isDynamic() ) {
      deallocate( begin_, capacity() );
   }
}
//*************************************************************************************************





//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::Reference
   SmallVector<T,N,A>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size(), "Invalid small vector access index" );
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference-to-const to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstReference
   SmallVector<T,N,A>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size(), "Invalid small vector access index" );
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid small vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::Reference
   SmallVector<T,N,A>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid small vector access index" );
   }
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid small vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstReference
   SmallVector<T,N,A>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid small vector access index" );
   }
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::Pointer
   SmallVector<T,N,A>::data() noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstPointer
   SmallVector<T,N,A>::data() const noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the small vector.
//
// \return Iterator to the first element of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::Iterator
   SmallVector<T,N,A>::begin() noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the small vector.
//
// \return Iterator to the first element of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstIterator
   SmallVector<T,N,A>::begin() const noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the small vector.
//
// \return Iterator to the first element of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstIterator
   SmallVector<T,N,A>::cbegin() const noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the small vector.
//
// \return Iterator just past the last element of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::Iterator
   SmallVector<T,N,A>::end() noexcept
{
   return end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the small vector.
//
// \return Iterator just past the last element of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstIterator
   SmallVector<T,N,A>::end() const noexcept
{
   return end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the small vector.
//
// \return Iterator just past the last element of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallVector<T,N,A>::ConstIterator
   SmallVector<T,N,A>::cend() const noexcept
{
   return end_;
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief List assignment to all vector elements.
//
// \param list The initializer list.
//
// This assignment operator offers the option to directly assign to all elements of the small
// vector by means of an initializer list:

   \code
   blaze::SmallVector<double,8UL> v;
   v = { 4.2, 6.3, -1.2 };
   \endcode

// The vector is resized according to the size of the initializer list and all its elements are
// assigned the values from the given initializer list.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
template< typename U >  // Type of the initializer list elements
inline SmallVector<T,N,A>& SmallVector<T,N,A>::operator=( initializer_list<U> list )
{
   resize( list.size() );
   std::copy( list.begin(), list.end(), begin_ );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for SmallVector.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
//
// The vector is resized according to the given small vector and initialized as a copy of this
// vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>& SmallVector<T,N,A>::operator=( const SmallVector& rhs )
{
   resize( rhs.size() );
   std::copy( rhs.begin(), rhs.end(), begin_ );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for SmallVector.
//
// \param rhs The vector to be moved into this instance.
// \return Reference to the assigned vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallVector<T,N,A>& SmallVector<T,N,A>::operator=( SmallVector&& rhs )
{
   resize( rhs.size() );
   std::move( rhs.begin(), rhs.end(), begin_ );

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the vector is empty.
//
// \return \a true in case the vector is empty, \a false if not.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline bool SmallVector<T,N,A>::empty() const noexcept
{
   return begin_ == end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size/dimension of the small vector.
//
// \return The current size of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline size_t SmallVector<T,N,A>::size() const noexcept
{
   return end_ - begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the small vector.
//
// \return The capacity of the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline size_t SmallVector<T,N,A>::capacity() const noexcept
{
   return final_ - begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the vector.
//
// \return void
//
// After the clear() function, the size of the vector is 0.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline void SmallVector<T,N,A>::clear()
{
   using blaze::destroy;

   destroy( begin_, end_ );

   if( isDynamic() ) {
      deallocate( begin_, capacity() );
   }

   begin_ = reinterpret_cast<T*>( v_ );
   end_   = begin_;
   final_ = begin_ + N;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
//
// \param n The new size of the vector.
// \return void
//
// This function resizes the vector using the given size to \a n. During this operation, new
// dynamic memory may be allocated in case the capacity of the vector is too small. New vector
// elements are not initialized!
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::resize( size_t n )
{
   using blaze::destroy;

   if( n > size() )
   {
      reserve( n );
      uninitialized_default_construct( begin_+size(), begin_+n );
      end_ = begin_ + n;
   }
   else if( n < size() )
   {
      destroy( begin_+size(), end_ );
      end_ = begin_ + n;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
//
// \param n The new size of the vector.
// \param value The initial value of new vector elements.
// \return void
//
// This function resizes the vector using the given size to \a n. During this operation, new
// dynamic memory may be allocated in case the capacity of the vector is too small. New vector
// elements are initialized to \a value!
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::resize( size_t n, const T& value )
{
   using blaze::destroy;

   if( n > size() )
   {
      reserve( n );
      std::uninitialized_fill( begin_+size(), begin_+n, value );
      end_ = begin_ + n;
   }
   else if( n < size() )
   {
      destroy( begin_+size(), end_ );
      end_ = begin_ + n;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the vector.
//
// \param n The new minimum capacity of the vector.
// \return void
//
// This function increases the capacity of the vector to at least \a n elements. The current
// values of the vector elements are preserved.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::reserve( size_t n )
{
   using blaze::destroy;

   const size_t oldCapacity( capacity() );

   if( n > oldCapacity )
   {
      const size_t oldSize( size() );
      T* tmp( allocate( n ) );

      if( IsNothrowMoveConstructible_v<T> ) {
         uninitialized_move( begin_, end_, tmp );
      }
      else {
         std::uninitialized_copy( begin_, end_, tmp );
      }

      destroy( begin_, end_ );

      if( isDynamic() ) {
         deallocate( begin_, oldCapacity );
      }

      final_ = tmp + n;
      end_   = tmp + oldSize;
      begin_ = tmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the vector by removing unused capacity. Please note
// that in case a reallocation occurs, all iterators (including end() iterators), all pointers
// and references to elements of this vector are invalidated.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::shrinkToFit()
{
   using blaze::destroy;

   const size_t oldCapacity( capacity() );
   const size_t oldSize    ( size() );

   if( isDynamic() && oldCapacity > oldSize )
   {
      T* tmp( allocate( oldSize ) );

      if( IsNothrowMoveConstructible_v<T> ) {
         uninitialized_move( begin_, end_, tmp );
      }
      else {
         std::uninitialized_copy( begin_, end_, tmp );
      }

      destroy( begin_, end_ );
      deallocate( begin_, oldCapacity );

      final_ = tmp + oldSize;
      end_   = final_;
      begin_ = tmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Adding an element to the end of the small vector.
//
// \param value The element to be added to the end of the small vector.
// \return void
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::pushBack( const T& value )
{
   using blaze::max;

   const size_t oldCapacity( capacity() );

   if( size() == oldCapacity ) {
      reserve( max( 2UL*oldCapacity, 7UL ) );
   }

   ::new ( end_ ) T( value );
   ++end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Adding an element to the end of the small vector.
//
// \param value The element to be added to the end of the small vector.
// \return void
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::pushBack( T&& value )
{
   using blaze::max;

   const size_t oldCapacity( capacity() );

   if( size() == oldCapacity ) {
      reserve( max( 2UL*oldCapacity, 7UL ) );
   }

   ::new ( end_ ) T( std::move( value ) );
   ++end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element at the specified position into the small vector.
//
// \param pos The position of the new element.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallVector<T,N,A>::Iterator
   SmallVector<T,N,A>::insert( Iterator pos, const T& value )
{
   using blaze::destroy;

   const size_t oldCapacity( capacity() );
   const size_t oldSize    ( size() );

   if( oldSize == oldCapacity )
   {
      const size_t newCapacity( 2UL*oldCapacity );
      const size_t index( pos - begin_ );

      T* tmp   ( allocate( newCapacity ) );
      T* newpos( tmp + index );

      uninitialized_move( begin_, pos, tmp );
      ::new ( newpos ) T( value );
      uninitialized_move( pos, end_, tmp+index+1UL );
      destroy( begin_, end_ );

      if( isDynamic() ) {
         deallocate( begin_, capacity() );
      }

      final_ = tmp + newCapacity;
      end_   = tmp + oldSize + 1UL;
      begin_ = tmp;

      return newpos;
   }
   else if( pos == end_ )
   {
      ::new( pos ) T( value );
      ++end_;
      return pos;
   }
   else
   {
      const auto tmp( end_ - 1UL );
      ::new ( end_ ) T( std::move( *tmp ) );

      try {
         std::move_backward( pos, tmp, end_ );
         destroy_at( pos );
         ::new ( pos ) T( value );
         ++end_;
      }
      catch( ... ) {
         destroy_at( end_ );
         throw;
      }

      return pos;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element at the specified position into the small vector.
//
// \param pos The position of the new element.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallVector<T,N,A>::Iterator
   SmallVector<T,N,A>::insert( Iterator pos, T&& value )
{
   using blaze::destroy;

   const size_t oldCapacity( capacity() );
   const size_t oldSize    ( size() );

   if( oldSize == oldCapacity )
   {
      const size_t newCapacity( 2UL*oldCapacity );
      const size_t index( pos - begin_ );

      T* tmp   ( allocate( newCapacity ) );
      T* newpos( tmp + index );

      uninitialized_move( begin_, pos, tmp );
      ::new ( newpos ) T( std::move( value ) );
      uninitialized_move( pos, end_, tmp+index+1UL );
      destroy( begin_, end_ );

      if( isDynamic() ) {
         deallocate( begin_, capacity() );
      }

      final_ = tmp + newCapacity;
      end_   = tmp + oldSize + 1UL;
      begin_ = tmp;

      return newpos;
   }
   else if( pos == end_ )
   {
      ::new( pos ) T( std::move( value ) );
      ++end_;
      return pos;
   }
   else
   {
      const auto tmp( end_ - 1UL );
      ::new ( end_ ) T( std::move( *tmp ) );

      try {
         std::move_backward( pos, tmp, end_ );
         destroy_at( pos );
         ::new ( pos ) T( std::move( value ) );
         ++end_;
      }
      catch( ... ) {
         destroy_at( end_ );
         throw;
      }

      return pos;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the small vector.
//
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallVector<T,N,A>::Iterator
   SmallVector<T,N,A>::erase( Iterator pos )
{
   std::move( pos+1UL, end_, pos );
   --end_;
   destroy_at( end_ );

   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the small vector.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the small vector.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallVector<T,N,A>::Iterator
   SmallVector<T,N,A>::erase( Iterator first, Iterator last )
{
   using blaze::destroy;

   BLAZE_USER_ASSERT( first <= last, "Invalid range detected" );

   const size_t n( last - first );

   std::move( last, end_, first );
   end_ -= n;
   destroy( end_, end_+n );

   return first;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two small vectors.
//
// \param sv The small vector to be swapped.
// \return void
//
// This function swaps the contents of two small vectors. Please note that this function is only
// guaranteed to not throw an exception if the move constructor of the underlying data type \a T
// is guaranteed to be noexcept.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallVector<T,N,A>::swap( SmallVector& sv ) noexcept( IsNothrowMoveConstructible_v<T> )
{
   using std::swap;
   using blaze::destroy;

   if( isDynamic() && sv.isDynamic() )
   {
      swap( begin_, sv.begin_ );
      swap( end_  , sv.end_   );
      swap( final_, sv.final_ );
   }
   else if( isDynamic() )
   {
      const size_t n( sv.size() );

      uninitialized_move( sv.begin_, sv.end_, reinterpret_cast<T*>( v_ ) );
      destroy( sv.begin_, sv.end_ );

      sv.begin_ = begin_;
      sv.end_   = end_;
      sv.final_ = final_;

      begin_ = reinterpret_cast<T*>( v_ );
      end_   = begin_ + n;
      final_ = begin_ + N;
   }
   else if( sv.isDynamic() )
   {
      const size_t n( size() );

      uninitialized_move( begin_, end_, reinterpret_cast<T*>( sv.v_ ) );
      destroy( begin_, end_ );

      begin_ = sv.begin_;
      end_   = sv.end_;
      final_ = sv.final_;

      sv.begin_ = reinterpret_cast<T*>( sv.v_ );
      sv.end_   = sv.begin_ + n;
      sv.final_ = sv.begin_ + N;
   }
   else if( size() > sv.size() )
   {
      const size_t n( size() - sv.size() );
      const auto pos = std::swap_ranges( sv.begin_, sv.end_, begin_ );

      uninitialized_move( pos, end_, sv.end_ );
      destroy( pos, end_ );
      end_    -= n;
      sv.end_ += n;
   }
   else
   {
      const size_t n( sv.size() - size() );
      const auto pos = std::swap_ranges( begin_, end_, sv.begin_ );

      uninitialized_move( pos, sv.end_, end_ );
      destroy( pos, sv.end_ );
      end_    += n;
      sv.end_ -= n;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the small vector uses its dynamic storage.
//
// \return \a true in case the dynamic storage is in use, \a false if not.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline bool SmallVector<T,N,A>::isDynamic() const noexcept
{
   return begin_ != reinterpret_cast<const T*>( v_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  SMALLVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SmallVector operators */
//@{
template< typename T1, size_t N1, typename A1, typename T2, size_t N2, typename A2 >
inline bool operator==( const SmallVector<T1,N1,A1>& lhs, const SmallVector<T2,N2,A2>& rhs );

template< typename T1, size_t N1, typename A1, typename T2, size_t N2, typename A2 >
inline bool operator!=( const SmallVector<T1,N1,A1>& lhs, const SmallVector<T2,N2,A2>& rhs );

template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::Iterator begin( SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator begin( const SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator cbegin( const SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::Iterator end( SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator end( const SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator cend( const SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline void clear( SmallVector<T,N,A>& sv );

template< typename T, size_t N, typename A >
inline void swap( SmallVector<T,N,A>& a, SmallVector<T,N,A>& b )
   noexcept( IsNothrowMoveConstructible_v<T> );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense vectors.
// \ingroup util
//
// \param lhs The left-hand side small vector for the comparison.
// \param rhs The right-hand side small vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
*/
template< typename T1    // Data type of the elements of the left-hand side vector
        , size_t N1      // Number of elements of the left-hand side vector
        , typename A1    // Type of the allocator of the left-hand side vector
        , typename T2    // Data type of the elements of the right-hand side vector
        , size_t N2      // Number of elements of the right-hand side vector
        , typename A2 >  // Type of the allocator of the right-hand side vector
inline bool operator==( const SmallVector<T1,N1,A1>& lhs, const SmallVector<T2,N2,A2>& rhs )
{
   if( lhs.size() != rhs.size() ) return false;

   return std::equal( lhs.begin(), lhs.end(), rhs.begin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense vectors.
// \ingroup util
//
// \param lhs The left-hand side small vector for the comparison.
// \param rhs The right-hand side small vector for the comparison.
// \return \a true if the two vectors are not equal, \a false if they are equal.
*/
template< typename T1    // Data type of the elements of the left-hand side vector
        , size_t N1      // Number of elements of the left-hand side vector
        , typename A1    // Type of the allocator of the left-hand side vector
        , typename T2    // Data type of the elements of the right-hand side vector
        , size_t N2      // Number of elements of the right-hand side vector
        , typename A2 >  // Type of the allocator of the right-hand side vector
inline bool operator!=( const SmallVector<T1,N1,A1>& lhs, const SmallVector<T2,N2,A2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given small vector.
// \ingroup util
//
// \param sv The given small vector.
// \return Iterator to the first element of the given small vector.
*/
template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::Iterator begin( SmallVector<T,N,A>& sv )
{
   return sv.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given small vector.
// \ingroup util
//
// \param sv The given small vector.
// \return Iterator to the first element of the given small vector.
*/
template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator begin( const SmallVector<T,N,A>& sv )
{
   return sv.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given small vector.
// \ingroup util
//
// \param sv The given small vector.
// \return Iterator to the first element of the given small vector.
*/
template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator cbegin( const SmallVector<T,N,A>& sv )
{
   return sv.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given small vector.
// \ingroup util
//
// \param sv The given small vector.
// \return Iterator just past the last element of the given small vector.
*/
template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::Iterator end( SmallVector<T,N,A>& sv )
{
   return sv.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given small vector.
// \ingroup util
//
// \param sv The given small vector.
// \return Iterator just past the last element of the given small vector.
*/
template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator end( const SmallVector<T,N,A>& sv )
{
   return sv.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given small vector.
// \ingroup util
//
// \param sv The given small vector.
// \return Iterator just past the last element of the given small vector.
*/
template< typename T, size_t N, typename A >
inline typename SmallVector<T,N,A>::ConstIterator cend( SmallVector<T,N,A>& sv )
{
   return sv.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given small vector.
// \ingroup util
//
// \param sv The small vector to be cleared.
// \return void
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline void clear( SmallVector<T,N,A>& sv )
{
   sv.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two small vectors.
// \ingroup util
//
// \param a The first vector to be swapped.
// \param b The second vector to be swapped.
// \return void
//
// This function swaps the contents of two small vectors. Please note that this function is only
// guaranteed to not throw an exception if the move constructor of the underlying data type \a T
// is guaranteed to be noexcept.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline void swap( SmallVector<T,N,A>& a, SmallVector<T,N,A>& b )
   noexcept( IsNothrowMoveConstructible_v<T> )
{
   a.swap( b );
}
//*************************************************************************************************

} // namespace blaze

#endif
