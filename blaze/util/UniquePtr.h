//=================================================================================================
/*!
//  \file blaze/util/UniquePtr.h
//  \brief Header file for the UniquePtr smart pointer class
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

#ifndef _BLAZE_UTIL_UNIQUEPTR_H_
#define _BLAZE_UTIL_UNIQUEPTR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Array.h>
#include <blaze/util/NonCopyable.h>
#include <blaze/util/Null.h>
#include <blaze/util/policies/PtrDelete.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scope-limited management of dynamically allocated resourses.
// \ingroup util
//
// The UniquePtr class implements a scope-restricted, lightweight smart pointer that manages a
// dynamically allocated resource. In contrast to other smart pointer implementations, UniquePtr
// is non-copyable and therefore restricted to manage resources within a single scope, but does
// so with a minimum of runtime overhead. The following example demonstrates the application of
// UniquePtr:

   \code
   {
      blaze::UniquePtr<std::string> mystring( new std::string( "My string" ) );

      // ... Working with the string

   } // The string is automatically destroyed at the end of scope according to the RAII principle
   \endcode

// Note that UniquePtr's interface is optimized for single objects created via the new operator
// and that by default DefaultDelete is used. For the management of dynamically allocated arrays,
// UniqueArray can be used:

   \code
   {
      blaze::UniqueArray<int> mystring( new int[10] );

      // ... Working with the integer array

   } // The array is automatically destroyed at the end of scope according to the RAII principle
   \endcode

// Note that the template argument \a T must not have any array extent (similar to the unique_ptr
// implementation of the C++11 standard, where the array extent differentiates between single
// objects and arrays). Providing a type with array extent will result in a compile time error!
*/
template< typename T                // Type of the resource
        , typename D = PtrDelete >  // Type of the deleter
class UniquePtr : private NonCopyable
{
 public:
   //**Type definitions****************************************************************************
   typedef typename RemoveReference<T>::Type*  Pointer;    //!< Pointer type of the managed resource.
   typedef typename RemoveReference<T>::Type&  Reference;  //!< Reference type of the managed resource.
   typedef D                                   Deleter;    //!< Type of the resource deleter.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline UniquePtr( Pointer ptr = NULL );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~UniquePtr();
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline Reference operator* () const /* throw() */;
   inline Pointer   operator->() const /* throw() */;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Pointer get() const /* throw() */;
   inline Pointer release() /* throw() */;
   inline void    reset( Pointer ptr = NULL ) /* throw() */;
   inline void    swap ( UniquePtr& up ) /* throw() */;
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
   Pointer ptr_;      //!< Pointer to the managed resource.
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
/*!\brief The default constructor for UniquePtr.
//
// \param ptr The resource to be managed by the unique pointer
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline UniquePtr<T,D>::UniquePtr( Pointer ptr )
   : ptr_    ( ptr       )  // Pointer to the managed resource
   , deleter_( Deleter() )  // Resource deleter.
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for UniquePtr.
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline UniquePtr<T,D>::~UniquePtr()
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
/*!\brief Direct access to the managed resource.
//
// \return Reference to the managed resource.
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline typename UniquePtr<T,D>::Reference UniquePtr<T,D>::operator*() const /* throw() */
{
   BLAZE_USER_ASSERT( ptr_, "Uninitialized unique pointer" );
   return *ptr_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the managed resource.
//
// \return Pointer to the managed resource.
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline typename UniquePtr<T,D>::Pointer UniquePtr<T,D>::operator->() const /* throw() */
{
   BLAZE_USER_ASSERT( ptr_, "Uninitialized unique pointer" );
   return ptr_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a pointer to the managed resource.
//
// \return Pointer to the managed resource or NULL if no resource is managed.
//
// This function returns a pointer to the managed resource (or NULL in case no resource is
// currently managed). Note however that the ownership remains with the unqiue pointer.
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline typename UniquePtr<T,D>::Pointer UniquePtr<T,D>::get() const /* throw() */
{
   return ptr_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Releases the ownership of the managed resource to the caller.
//
// \return Pointer to the managed resource or NULL if no resource is managed.
//
// This function returns a pointer to the managed resource (or NULL in case no resource is
// currently managed). The ownership of the resource is released and passed to the caller.
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline typename UniquePtr<T,D>::Pointer UniquePtr<T,D>::release() /* throw() */
{
   Pointer tmp( ptr_ );
   ptr_ = NULL;
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resets the unique pointer and replaces the managed resource with the given resource.
//
// \param ptr The new resource to be managed by the unique pointer.
// \return void
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline void UniquePtr<T,D>::reset( Pointer ptr ) /* throw() */
{
   if( ptr != ptr_ ) {
      UniquePtr( ptr ).swap( *this );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two unique pointers.
//
// \param ptr The unique pointer to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline void UniquePtr<T,D>::swap( UniquePtr& ptr ) /* throw() */
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
template< typename T    // Type of the resource
        , typename D >  // Type of the deleter
inline UniquePtr<T,D>::operator bool() const /* throw() */
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
/*!\name UniquePtr operators */
//@{
template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator==( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator!=( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator<( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator<=( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator>( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs );

template< typename T1, typename D1, typename T2, typename D2 >
inline bool operator>=( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs );

template< typename T, typename D >
inline bool operator==( const UniquePtr<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator!=( const UniquePtr<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator<( const UniquePtr<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator>( const UniquePtr<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator<=( const UniquePtr<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator>=( const UniquePtr<T,D>& ptr, const Null& null );

template< typename T, typename D >
inline bool operator==( const Null& null, const UniquePtr<T,D>& ptr );

template< typename T, typename D >
inline bool operator!=( const Null& null, const UniquePtr<T,D>& ptr );

template< typename T, typename D >
inline bool operator<( const Null& null, const UniquePtr<T,D>& ptr );

template< typename T, typename D >
inline bool operator>( const Null& null, const UniquePtr<T,D>& ptr );

template< typename T, typename D >
inline bool operator<=( const Null& null, const UniquePtr<T,D>& ptr );

template< typename T, typename D >
inline bool operator>=( const Null& null, const UniquePtr<T,D>& ptr );

template< typename T, typename D >
inline void swap( UniquePtr<T,D>& a, UniquePtr<T,D>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two UniquePtr objects.
//
// \param lhs The left-hand side unique pointer.
// \param rhs The right-hand side unique pointer.
// \return \a true if the two pointers refer to the same element, \a false if they don't.
*/
template< typename T1    // Resource type of the left-hand side unique pointer
        , typename D1    // Deleter type of the left-hand side unique pointer
        , typename T2    // Resource type of the right-hand side unique pointer
        , typename D2 >  // Deleter type of the right-hand side unique pointer
inline bool operator==( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs )
{
   return lhs.get() == rhs.get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two UniquePtr objects.
//
// \param lhs The left-hand side unique pointer.
// \param rhs The right-hand side unique pointer.
// \return \a true if the two pointers do not refer to the same element, \a false if they do.
*/
template< typename T1    // Resource type of the left-hand side unique pointer
        , typename D1    // Deleter type of the left-hand side unique pointer
        , typename T2    // Resource type of the right-hand side unique pointer
        , typename D2 >  // Deleter type of the right-hand side unique pointer
inline bool operator!=( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs )
{
   return lhs.get() != rhs.get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two UniquePtr objects.
//
// \param lhs The left-hand side unique pointer.
// \param rhs The right-hand side unique pointer.
// \return \a true if the left pointer is less than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique pointer
        , typename D1    // Deleter type of the left-hand side unique pointer
        , typename T2    // Resource type of the right-hand side unique pointer
        , typename D2 >  // Deleter type of the right-hand side unique pointer
inline bool operator<( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs )
{
   return lhs.get() < rhs.get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two UniquePtr objects.
//
// \param lhs The left-hand side unique pointer.
// \param rhs The right-hand side unique pointer.
// \return \a true if the left pointer is greater than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique pointer
        , typename D1    // Deleter type of the left-hand side unique pointer
        , typename T2    // Resource type of the right-hand side unique pointer
        , typename D2 >  // Deleter type of the right-hand side unique pointer
inline bool operator>( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs )
{
   return rhs < lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two UniquePtr objects.
//
// \param lhs The left-hand side unique pointer.
// \param rhs The right-hand side unique pointer.
// \return \a true if the left pointer is less or equal than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique pointer
        , typename D1    // Deleter type of the left-hand side unique pointer
        , typename T2    // Resource type of the right-hand side unique pointer
        , typename D2 >  // Deleter type of the right-hand side unique pointer
inline bool operator<=( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs )
{
   return !( rhs < lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two UniquePtr objects.
//
// \param lhs The left-hand side unique pointer.
// \param rhs The right-hand side unique pointer.
// \return \a true if the left pointer is greater or equal than the right pointer, \a false if not.
*/
template< typename T1    // Resource type of the left-hand side unique pointer
        , typename D1    // Deleter type of the left-hand side unique pointer
        , typename T2    // Resource type of the right-hand side unique pointer
        , typename D2 >  // Deleter type of the right-hand side unique pointer
inline bool operator>=( const UniquePtr<T1,D1>& lhs, const UniquePtr<T2,D2>& rhs )
{
   return !( lhs < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a UniquePtr and NULL.
//
// \param ptr The left-hand side unique pointer.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique pointer is NULL, \a false if it isn't.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator==( const UniquePtr<T,D>& ptr, const Null& null )
{
   return ptr.get() == null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a UniquePtr and NULL.
//
// \param ptr The left-hand side unique pointer.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique pointer is not NULL, \a false if it is.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator!=( const UniquePtr<T,D>& ptr, const Null& null )
{
   return !( ptr == null );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a UniquePtr and NULL.
//
// \param ptr The left-hand side unique pointer.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique pointer is less than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator<( const UniquePtr<T,D>& ptr, const Null& null )
{
   return ptr.get() < null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a UniquePtr and NULL.
//
// \param ptr The left-hand side unique pointer.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique pointer is greater than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator>( const UniquePtr<T,D>& ptr, const Null& null )
{
   return ptr.get() > null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a UniquePtr and NULL.
//
// \param ptr The left-hand side unique pointer.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique pointer is less or equal than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator<=( const UniquePtr<T,D>& ptr, const Null& null )
{
   return !( ptr > null );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a UniquePtr and NULL.
//
// \param ptr The left-hand side unique pointer.
// \param null The right-hand side NULL pointer.
// \return \a true if the unique pointer is greater or equal than NULL, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator>=( const UniquePtr<T,D>& ptr, const Null& null )
{
   return !( ptr < null );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between NULL and a UniquePtr.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique pointer.
// \return \a true if the unique pointer is NULL, \a false if it isn't.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator==( const Null& null, const UniquePtr<T,D>& ptr )
{
   return ptr == null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between NULL and a UniquePtr.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique pointer.
// \return \a true if the unique pointer is not NULL, \a false if it is.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator!=( const Null& null, const UniquePtr<T,D>& ptr )
{
   return ptr != null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between NULL and a UniquePtr.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique pointer.
// \return \a true if NULL is less than the unique pointer, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator<( const Null& null, const UniquePtr<T,D>& ptr )
{
   return ptr > null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between NULL and a UniquePtr.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique pointer.
// \return \a true if NULL is greater than the unique pointer, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator>( const Null& null, const UniquePtr<T,D>& ptr )
{
   return ptr < null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between NULL and a UniquePtr.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique pointer.
// \return \a true if NULL is less or equal than the unique pointer, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator<=( const Null& null, const UniquePtr<T,D>& ptr )
{
   return ptr >= null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between NULL and a UniquePtr.
//
// \param null The left-hand side NULL pointer.
// \param ptr The right-hand side unique pointer.
// \return \a true if NULL is greater or equal than the unique pointer, \a false if not.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline bool operator>=( const Null& null, const UniquePtr<T,D>& ptr )
{
   return ptr <= null;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two unique pointers.
//
// \param a The first unique pointer to be swapped.
// \param b The second unique pointer to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename T    // Resource type of the unique pointer
        , typename D >  // Deleter type of the unique pointer
inline void swap( UniquePtr<T,D>& a, UniquePtr<T,D>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************

} // namespace blaze

#endif
