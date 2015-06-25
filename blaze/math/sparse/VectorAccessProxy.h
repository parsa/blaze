//=================================================================================================
/*!
//  \file blaze/math/sparse/VectorAccessProxy.h
//  \brief Header file for the VectorAccessProxy class
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

#ifndef _BLAZE_MATH_SPARSE_VECTORACCESSPROXY_H_
#define _BLAZE_MATH_SPARSE_VECTORACCESSPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <ostream>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/proxy/Proxy.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for sparse, N-dimensional vectors.
// \ingroup math
//
// The VectorAccessProxy provides safe access to the elements of a non-const sparse vector.\n
// The proxied access to the elements of a sparse vector is necessary since it may be possible
// that several insertion operations happen in the same statement. The following code illustrates
// this with two examples by means of the CompressedVector class:

   \code
   blaze::CompressedVector<double> a( 5 );

   // Standard usage of the subscript operator to initialize a vector element.
   // Only a single vector matrix element is accessed!
   a[0] = 1.0;

   // Initialization of a vector element via another vector element.
   // Two sparse vector accesses in one statement!
   a[1] = a[0];

   // Multiple accesses to elements of the sparse vector in one statement!
   const double result = a[0] + a[2] + a[4];
   \endcode

// The problem (especially with the last statement) is that several insertion operations might
// take place due to the access via the subscript operator. If the subscript operator would
// return a direct reference to one of the accessed elements, this reference might be invalidated
// during the evaluation of a subsequent subscript operator, which results in undefined behavior.
// This class provides the necessary functionality to guarantee a safe access to the sparse vector
// elements while preserving the intuitive use of the subscript operator.
//
*/
template< typename VT >  // Type of the sparse vector
class VectorAccessProxy : public Proxy< VectorAccessProxy<VT>, typename VT::ElementType >
{
 public:
   //**Type definitions****************************************************************************
   typedef typename VT::ElementType  RepresentedType;  //!< Type of the represented sparse vector element.
   typedef RepresentedType&          RawReference;     //!< Raw reference to the represented element.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline VectorAccessProxy( VT& sv, size_t i );
            inline VectorAccessProxy( const VectorAccessProxy& vap );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~VectorAccessProxy();
   //@}
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                          inline const VectorAccessProxy& operator= ( const VectorAccessProxy& vap ) const;
   template< typename T > inline const VectorAccessProxy& operator= ( const T& value ) const;
   template< typename T > inline const VectorAccessProxy& operator+=( const T& value ) const;
   template< typename T > inline const VectorAccessProxy& operator-=( const T& value ) const;
   template< typename T > inline const VectorAccessProxy& operator*=( const T& value ) const;
   template< typename T > inline const VectorAccessProxy& operator/=( const T& value ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline RawReference get()          const;
   inline bool         isRestricted() const;
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator RawReference() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT&    sv_;  //!< Reference to the accessed sparse vector.
   size_t i_;   //!< Index of the accessed sparse vector element.
   //@}
   //**********************************************************************************************

   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   void* operator&() const;  //!< Address operator (private & undefined)
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
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
/*!\brief Initialization constructor for a VectorAccessProxy.
//
// \param sv Reference to the accessed sparse vector.
// \param i The index of the accessed sparse vector element.
*/
template< typename VT >  // Type of the sparse vector
inline VectorAccessProxy<VT>::VectorAccessProxy( VT& sv, size_t i )
   : sv_( sv )  // Reference to the accessed sparse vector
   , i_ ( i  )  // Index of the accessed sparse vector element
{
   const typename VT::Iterator element( sv_.find( i_ ) );
   if( element == sv_.end() )
      sv_.insert( i_, RepresentedType() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for VectorAccessProxy.
//
// \param vap Sparse vector access proxy to be copied.
*/
template< typename VT >  // Type of the sparse vector
inline VectorAccessProxy<VT>::VectorAccessProxy( const VectorAccessProxy& vap )
   : sv_( vap.sv_ )  // Reference to the accessed sparse vector
   , i_ ( vap.i_  )  // Index of the accessed sparse vector element
{
   BLAZE_INTERNAL_ASSERT( sv_.find( i_ ) != sv_.end(), "Missing vector element detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for VectorAccessProxy.
*/
template< typename VT >  // Type of the sparse vector
inline VectorAccessProxy<VT>::~VectorAccessProxy()
{
   const typename VT::Iterator element( sv_.find( i_ ) );
   if( element != sv_.end() && isDefault( element->value() ) )
      sv_.erase( element );
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for VectorAccessProxy.
//
// \param vap Sparse vector access proxy to be copied.
// \return Reference to the assigned access proxy.
*/
template< typename VT >  // Type of the sparse vector
inline const VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator=( const VectorAccessProxy& vap ) const
{
   get() = vap.get();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the accessed sparse vector element.
//
// \param value The new value of the sparse vector element.
// \return Reference to the assigned access proxy.
*/
template< typename VT >  // Type of the sparse vector
template< typename T >   // Type of the right-hand side value
inline const VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator=( const T& value ) const
{
   get() = value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the accessed sparse vector element.
//
// \param value The right-hand side value to be added to the sparse vector element.
// \return Reference to the assigned access proxy.
*/
template< typename VT >  // Type of the sparse vector
template< typename T >   // Type of the right-hand side value
inline const VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator+=( const T& value ) const
{
   get() += value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the accessed sparse vector element.
//
// \param value The right-hand side value to be subtracted from the sparse vector element.
// \return Reference to the assigned access proxy.
*/
template< typename VT >  // Type of the sparse vector
template< typename T >   // Type of the right-hand side value
inline const VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator-=( const T& value ) const
{
   get() -= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the accessed sparse vector element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned access proxy.
*/
template< typename VT >  // Type of the sparse vector
template< typename T >   // Type of the right-hand side value
inline const VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator*=( const T& value ) const
{
   get() *= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the accessed sparse vector element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned access proxy.
*/
template< typename VT >  // Type of the sparse vector
template< typename T >   // Type of the right-hand side value
inline const VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator/=( const T& value ) const
{
   get() /= value;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returning the value of the accessed sparse vector element.
//
// \return Direct/raw reference to the accessed sparse vector element.
*/
template< typename VT >  // Type of the sparse vector
inline typename VectorAccessProxy<VT>::RawReference VectorAccessProxy<VT>::get() const
{
   const typename VT::Iterator element( sv_.find( i_ ) );
   BLAZE_INTERNAL_ASSERT( element != sv_.end(), "Missing vector element detected" );
   return element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the proxy represents a restricted sparse vector element..
//
// \return \a true in case access to the sparse vector element is restricted, \a false if not.
*/
template< typename VT >  // Type of the sparse vector
inline bool VectorAccessProxy<VT>::isRestricted() const
{
   return false;
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion to the accessed sparse vector element.
//
// \return Direct/raw reference to the accessed sparse vector element.
*/
template< typename VT >  // Type of the sparse vector
inline VectorAccessProxy<VT>::operator RawReference() const
{
   return get();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name VectorAccessProxy operators */
//@{
template< typename VT1, typename VT2 >
inline bool operator==( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs );

template< typename VT, typename T >
inline bool operator==( const VectorAccessProxy<VT>& lhs, const T& rhs );

template< typename T, typename VT >
inline bool operator==( const T& lhs, const VectorAccessProxy<VT>& rhs );

template< typename VT1, typename VT2 >
inline bool operator!=( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs );

template< typename VT, typename T >
inline bool operator!=( const VectorAccessProxy<VT>& lhs, const T& rhs );

template< typename T, typename VT >
inline bool operator!=( const T& lhs, const VectorAccessProxy<VT>& rhs );

template< typename VT1, typename VT2 >
inline bool operator<( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs );

template< typename VT, typename T >
inline bool operator<( const VectorAccessProxy<VT>& lhs, const T& rhs );

template< typename T, typename VT >
inline bool operator<( const T& lhs, const VectorAccessProxy<VT>& rhs );

template< typename VT1, typename VT2 >
inline bool operator>( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs );

template< typename VT, typename T >
inline bool operator>( const VectorAccessProxy<VT>& lhs, const T& rhs );

template< typename T, typename VT >
inline bool operator>( const T& lhs, const VectorAccessProxy<VT>& rhs );

template< typename VT1, typename VT2 >
inline bool operator<=( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs );

template< typename VT, typename T >
inline bool operator<=( const VectorAccessProxy<VT>& lhs, const T& rhs );

template< typename T, typename VT >
inline bool operator<=( const T& lhs, const VectorAccessProxy<VT>& rhs );

template< typename VT1, typename VT2 >
inline bool operator>=( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs );

template< typename VT, typename T >
inline bool operator>=( const VectorAccessProxy<VT>& lhs, const T& rhs );

template< typename T, typename VT >
inline bool operator>=( const T& lhs, const VectorAccessProxy<VT>& rhs );

template< typename VT >
inline std::ostream& operator<<( std::ostream& os, const VectorAccessProxy<VT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two VectorAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename VT1, typename VT2 >
inline bool operator==( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs )
{
   return ( lhs.get() == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a VectorAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename VT, typename T >
inline bool operator==( const VectorAccessProxy<VT>& lhs, const T& rhs )
{
   return ( lhs.get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a VectorAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename VT >
inline bool operator==( const T& lhs, const VectorAccessProxy<VT>& rhs )
{
   return ( lhs == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two VectorAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename VT1, typename VT2 >
inline bool operator!=( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs )
{
   return ( lhs.get() != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a VectorAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename VT, typename T >
inline bool operator!=( const VectorAccessProxy<VT>& lhs, const T& rhs )
{
   return ( lhs.get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a VectorAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename VT >
inline bool operator!=( const T& lhs, const VectorAccessProxy<VT>& rhs )
{
   return ( lhs != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two VectorAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename VT1, typename VT2 >
inline bool operator<( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs )
{
   return ( lhs.get() < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a VectorAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename VT, typename T >
inline bool operator<( const VectorAccessProxy<VT>& lhs, const T& rhs )
{
   return ( lhs.get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a VectorAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename VT >
inline bool operator<( const T& lhs, const VectorAccessProxy<VT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two VectorAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename VT1, typename VT2 >
inline bool operator>( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs )
{
   return ( lhs.get() > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a VectorAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename VT, typename T >
inline bool operator>( const VectorAccessProxy<VT>& lhs, const T& rhs )
{
   return ( lhs.get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a VectorAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename VT >
inline bool operator>( const T& lhs, const VectorAccessProxy<VT>& rhs )
{
   return ( lhs > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two VectorAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename VT1, typename VT2 >
inline bool operator<=( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs )
{
   return ( lhs.get() <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a VectorAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename VT, typename T >
inline bool operator<=( const VectorAccessProxy<VT>& lhs, const T& rhs )
{
   return ( lhs.get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a VectorAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename VT >
inline bool operator<=( const T& lhs, const VectorAccessProxy<VT>& rhs )
{
   return ( lhs <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two VectorAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename VT1, typename VT2 >
inline bool operator>=( const VectorAccessProxy<VT1>& lhs, const VectorAccessProxy<VT2>& rhs )
{
   return ( lhs.get() >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a VectorAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side VectorAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename VT, typename T >
inline bool operator>=( const VectorAccessProxy<VT>& lhs, const T& rhs )
{
   return ( lhs.get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a VectorAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side VectorAccessProxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename VT >
inline bool operator>=( const T& lhs, const VectorAccessProxy<VT>& rhs )
{
   return ( lhs >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the VectorAccessProxy class template.
// \ingroup math
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant proxy object.
// \return Reference to the output stream.
*/
template< typename VT >
inline std::ostream& operator<<( std::ostream& os, const VectorAccessProxy<VT>& proxy )
{
   return os << proxy.get();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name VectorAccessProxy global functions */
//@{
template< typename VT >
inline void reset( const VectorAccessProxy<VT>& proxy );

template< typename VT >
inline void clear( const VectorAccessProxy<VT>& proxy );

template< typename VT >
inline bool isDefault( const VectorAccessProxy<VT>& proxy );

template< typename VT >
inline void swap( const VectorAccessProxy<VT>& a, const VectorAccessProxy<VT>& b ) /* throw() */;

template< typename VT, typename T >
inline void swap( const VectorAccessProxy<VT>& a, T& b ) /* throw() */;

template< typename T, typename VT >
inline void swap( T& a, const VectorAccessProxy<VT>& v ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the represented element to the default initial values.
// \ingroup math
//
// \param proxy The given access proxy.
// \return void
*/
template< typename VT >
inline void reset( const VectorAccessProxy<VT>& proxy )
{
   using blaze::reset;

   reset( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
// \ingroup math
//
// \param proxy The given access proxy.
// \return void
*/
template< typename VT >
inline void clear( const VectorAccessProxy<VT>& proxy )
{
   using blaze::clear;

   clear( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is in default state.
// \ingroup math
//
// \param proxy The given access proxy.
// \return \a true in case the represented element is in default state, \a false otherwise.
//
// This function checks whether the element represented by the access proxy is in default state.
// In case it is in default state, the function returns \a true, otherwise it returns \a false.
*/
template< typename VT >
inline bool isDefault( const VectorAccessProxy<VT>& proxy )
{
   using blaze::isDefault;

   return isDefault( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two access proxies.
// \ingroup math
//
// \param a The first access proxy to be swapped.
// \param b The second access proxy to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename VT >
inline void swap( const VectorAccessProxy<VT>& a, const VectorAccessProxy<VT>& b ) /* throw() */
{
   using std::swap;

   swap( a.get(), b.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of an access proxy with another element.
// \ingroup math
//
// \param a The access proxy to be swapped.
// \param b The other element to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename VT, typename T >
inline void swap( const VectorAccessProxy<VT>& a, T& b ) /* throw() */
{
   using std::swap;

   swap( a.get(), b );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of an access proxy with another element.
// \ingroup math
//
// \param a The other element to be swapped.
// \param b The access proxy to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename T, typename VT >
inline void swap( T& a, const VectorAccessProxy<VT>& b ) /* throw() */
{
   using std::swap;

   swap( a, b.get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
