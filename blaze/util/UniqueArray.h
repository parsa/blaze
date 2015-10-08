//=================================================================================================
/*!
//  \file blaze/util/UniqueArray.h
//  \brief Header file for the UniqueArray smart pointer class
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

#ifndef _BLAZE_UTIL_UNIQUEARRAY_H_
#define _BLAZE_UTIL_UNIQUEARRAY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Array.h>
#include <blaze/util/NonCopyable.h>
#include <blaze/util/Null.h>
#include <blaze/util/policies/ArrayDelete.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scope-limited management of dynamically allocated arrays.
// \ingroup util
//
// The UniqueArray class implements a scope-restricted, lightweight smart pointer that manages a
// dynamically allocated array. In contrast to other smart pointer implementations, UniqueArray
// is non-copyable and therefore restricted to manage arrays within a single scope, but does so
// so with a minimum of runtime overhead. The following example demonstrates the application of
// UniqueArray:

   \code
   {
      blaze::UniqueArray<int> mystring( new int[10] );

      // ... Working with the integer array

   } // The array is automatically destroyed at the end of scope according to the RAII principle
   \endcode

// Note that UniqueArray's interface is optimized for arrays and that by default ArrayDelete is
// used. For the management of single objects, UniquePtr can be used:

   \code
   {
      blaze::UniquePtr<std::string> mystring( new std::string( "My string" ) );

      // ... Working with the string

   } // The string is automatically destroyed at the end of scope according to the RAII principle
   \endcode

// Note that the template argument \a T must not have any array extent (similar to the unique_ptr
// implementation of the C++11 standard, where the array extent differentiates between single
// objects and arrays). Providing a type with array extent will result in a compile time error!
*/
template< typename T                  // Type of the array elements
        , typename D = ArrayDelete >  // Type of the deleter
class UniqueArray : private NonCopyable
{
 public:
   //**Type definitions****************************************************************************
   typedef typename RemoveReference<T>::Type*  Pointer;    //!< Pointer type of the managed array elements.
   typedef typename RemoveReference<T>::Type&  Reference;  //!< Reference type of the managed array elements.
   typedef D                                   Deleter;    //!< Type of the resource deleter.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline UniqueArray( Pointer ptr = NULL );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~UniqueArray();
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline Reference operator[]( size_t index ) const  /* throw() */;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Pointer get() const /* throw() */;
   inline Pointer release() /* throw() */;
   inline void    reset( Pointer ptr = NULL ) /* throw() */;
   inline void    swap ( UniqueArray& up ) /* throw() */;
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator bool() const /* throw() */;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Pointer ptr_;      //!< Pointer to the managed array.
   Deleter deleter_;  //!< Resource deleter.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_ARRAY_TYPE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for the UniqueArray specialization.
//
// \param ptr The array to be managed by the unique array
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline UniqueArray<T,D>::UniqueArray( Pointer ptr )
   : ptr_    ( ptr       )  // Pointer to the managed array
   , deleter_( Deleter() )  // Resource deleter
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for the UniqueArray specialization.
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline UniqueArray<T,D>::~UniqueArray()
{
   deleter_( ptr_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ACCESS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Direct access to the array elements.
//
// \return Reference to the managed array.
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline typename UniqueArray<T,D>::Reference UniqueArray<T,D>::operator[]( size_t index ) const /* throw() */
{
   BLAZE_USER_ASSERT( ptr_, "Uninitialized unique pointer" );
   return ptr_[index];
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a pointer to the managed array.
//
// \return Pointer to the managed array or NULL if no array is managed.
//
// This function returns a pointer to the managed array (or NULL in case no array is
// currently managed). Note however that the ownership remains with the unqiue pointer.
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline typename UniqueArray<T,D>::Pointer UniqueArray<T,D>::get() const /* throw() */
{
   return ptr_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Releases the ownership of the managed array to the caller.
//
// \return Pointer to the managed array or NULL if no array is managed.
//
// This function returns a pointer to the managed array (or NULL in case no array is
// currently managed). The ownership of the array is released and passed to the caller.
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline typename UniqueArray<T,D>::Pointer UniqueArray<T,D>::release() /* throw() */
{
   Pointer tmp( ptr_ );
   ptr_ = NULL;
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resets the unique array and replaces the managed array with the given array.
//
// \param ptr The new array to be managed by the unique array.
// \return void
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline void UniqueArray<T,D>::reset( Pointer ptr ) /* throw() */
{
   if( ptr != ptr_ ) {
      UniqueArray( ptr ).swap( *this );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two unique arrays.
//
// \param ptr The unique array to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline void UniqueArray<T,D>::swap( UniqueArray& ptr ) /* throw() */
{
   Pointer tmp( ptr_ );
   ptr_ = ptr.ptr_;
   ptr.ptr_ = tmp;
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the unique pointer is set to a non-NULL pointer.
//
// \return \a true in case the unique pointer is not NULL, \a false if it is NULL.
*/
template< typename T    // Type of the array elements
        , typename D >  // Type of the deleter
inline UniqueArray<T,D>::operator bool() const /* throw() */
{
   return ( ptr_ != NULL );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name UniqueArray operators */
//@{
template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator==( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator!=( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator<( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator<=( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator>( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator>=( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs );

template< typename T, typename D >
inline bool operator==( const UniqueArray<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator!=( const UniqueArray<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator<( const UniqueArray<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator>( const UniqueArray<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator<=( const UniqueArray<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator>=( const UniqueArray<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator==( const Null& null, const UniqueArray<T,D>& ptr );

template< typename T, typename D >
inline bool operator!=( const Null& null, const UniqueArray<T,D>& ptr );

template< typename T, typename D >
inline bool operator<( const Null& null, const UniqueArray<T,D>& ptr );

template< typename T, typename D >
inline bool operator>( const Null& null, const UniqueArray<T,D>& ptr );

template< typename T, typename D >
inline bool operator<=( const Null& null, const UniqueArray<T,D>& ptr );

template< typename T, typename D >
inline bool operator>=( const Null& null, const UniqueArray<T,D>& ptr );

template< typename T, typename D >
inline void swap( UniqueArray<T,D>& a, UniqueArray<T,D>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two UniqueArray objects.
//
// \param lhs The left-hand side unique array.
// \param rhs The right-hand side unique array.
// \return \a true if the two pointers refer to the same element, \a false if they don't.
*/
template< typename T1    // Resource type of the left-hand side unique array
        , typename D1    // Deleter type of the left-hand side unique array
        , typename T2    // Resource type of the right-hand side unique array
        , typename D2 >  // Deleter type of the right-hand side unique array
inline bool operator==( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs )
{
   return lhs.get() == rhs.get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two UniqueArray objects.
//
// \param lhs The left-hand side unique array.
// \param rhs The right-hand side unique array.
// \return \a true if the two pointers do not refer to the same element, \a false if they do.
*/
template< typename T1    // Resource type of the left-hand side unique array
        , typename D1    // Deleter type of the left-hand side unique array
        , typename T2    // Resource type of the right-hand side unique array
        , typename D2 >  // Deleter type of the right-hand side unique array
inline bool operator!=( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs )
{
   return lhs.get() != rhs.get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two UniqueArray objects.
//
// \param lhs The left-hand side unique array.
// \param rhs The right-hand side unique array.
// \return \a true if the left pointer is less than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique array
        , typename D1    // Deleter type of the left-hand side unique array
        , typename T2    // Resource type of the right-hand side unique array
        , typename D2 >  // Deleter type of the right-hand side unique array
inline bool operator<( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs )
{
   return lhs.get() < rhs.get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two UniqueArray objects.
//
// \param lhs The left-hand side unique array.
// \param rhs The right-hand side unique array.
// \return \a true if the left pointer is greater than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique array
        , typename D1    // Deleter type of the left-hand side unique array
        , typename T2    // Resource type of the right-hand side unique array
        , typename D2 >  // Deleter type of the right-hand side unique array
inline bool operator>( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs )
{
   return rhs < lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two UniqueArray objects.
//
// \param lhs The left-hand side unique array.
// \param rhs The right-hand side unique array.
// \return \a true if the left pointer is less or equal than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique array
        , typename D1    // Deleter type of the left-hand side unique array
        , typename T2    // Resource type of the right-hand side unique array
        , typename D2 >  // Deleter type of the right-hand side unique array
inline bool operator<=( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs )
{
   return !( rhs < lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two UniqueArray objects.
//
// \param lhs The left-hand side unique array.
// \param rhs The right-hand side unique array.
// \return \a true if the left pointer is greater or equal than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique array
        , typename D1    // Deleter type of the left-hand side unique array
        , typename T2    // Resource type of the right-hand side unique array
        , typename D2 >  // Deleter type of the right-hand side unique array
inline bool operator>=( const UniqueArray<T1,D1>& lhs, const UniqueArray<T2,D2>& rhs )
{
   return !( lhs < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a UniqueArray and NULL.
//
// \param ptr The left-hand side unique array.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique array is NULL, \a false if it isn't.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator==( const UniqueArray<T,D>& ptr, const Null& null )
{
   return ptr.get() == null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a UniqueArray and NULL.
//
// \param ptr The left-hand side unique array.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique array is not NULL, \a false if it is.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator!=( const UniqueArray<T,D>& ptr, const Null& null )
{
   return !( ptr == null );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a UniqueArray and NULL.
//
// \param ptr The left-hand side unique array.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique array is less than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator<( const UniqueArray<T,D>& ptr, const Null& null )
{
   return ptr.get() < null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a UniqueArray and NULL.
//
// \param ptr The left-hand side unique array.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique array is greater than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator>( const UniqueArray<T,D>& ptr, const Null& null )
{
   return ptr.get() > null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a UniqueArray and NULL.
//
// \param ptr The left-hand side unique array.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique array is less or equal than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator<=( const UniqueArray<T,D>& ptr, const Null& null )
{
   return !( ptr > null );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a UniqueArray and NULL.
//
// \param ptr The left-hand side unique array.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique array is greater or equal than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator>=( const UniqueArray<T,D>& ptr, const Null& null )
{
   return !( ptr < null );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between NULL and a UniqueArray.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique array.
// \return \a true if the unique array is NULL, \a false if it isn't.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator==( const Null& null, const UniqueArray<T,D>& ptr )
{
   return ptr == null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between NULL and a UniqueArray.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique array.
// \return \a true if the unique array is not NULL, \a false if it is.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator!=( const Null& null, const UniqueArray<T,D>& ptr )
{
   return ptr != null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between NULL and a UniqueArray.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique array.
// \return \a true if NULL is less than the unique array, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator<( const Null& null, const UniqueArray<T,D>& ptr )
{
   return ptr > null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between NULL and a UniqueArray.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique array.
// \return \a true if NULL is greater than the unique array, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator>( const Null& null, const UniqueArray<T,D>& ptr )
{
   return ptr < null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between NULL and a UniqueArray.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique array.
// \return \a true if NULL is less or equal than the unique array, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator<=( const Null& null, const UniqueArray<T,D>& ptr )
{
   return ptr >= null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between NULL and a UniqueArray.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique array.
// \return \a true if NULL is greater or equal than the unique array, \a false if not.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline bool operator>=( const Null& null, const UniqueArray<T,D>& ptr )
{
   return ptr <= null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two unique arrays.
//
// \param a The first unique array to be swapped.
// \param b The second unique array to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename T    // Resource type of the unique array
        , typename D >  // Deleter type of the unique array
inline void swap( UniqueArray<T,D>& a, UniqueArray<T,D>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************

} // namespace blaze

#endif
