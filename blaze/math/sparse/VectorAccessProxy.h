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

#include <blaze/math/shims/IsDefault.h>
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
   CompressedVector<real> a( 5 );

   // Standard usage of the subscript operator to initialize a vector element.
   // Only a single vector matrix element is accessed!
   a[0] = 1.0;

   // Initialization of a vector element via another vector element.
   // Two sparse vector accesses in one statement!
   a[1] = a[0];

   // Multiple accesses to elements of the sparse vector in one statement!
   const real result = a[0] + a[2] + a[4];
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
class VectorAccessProxy
{
 public:
   //**Type definitions****************************************************************************
   typedef VT                        VectorType;      //!< Type of the accessed sparse vector.
   typedef typename VT::ElementType  ElementType;     //!< Type of the accessed sparse vector element.
   typedef ElementType&              Reference;       //!< Reference type of the accessed element.
   typedef const ElementType&        ConstReference;  //!< Reference type of the accessed constant element.
   typedef typename VT::Iterator     Iterator;        //!< Iterator type of the accessed sparse vector.
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
                          inline VectorAccessProxy& operator= ( const VectorAccessProxy& vap );
   template< typename T > inline VectorAccessProxy& operator= ( const T& value );
   template< typename T > inline VectorAccessProxy& operator+=( const T& value );
   template< typename T > inline VectorAccessProxy& operator-=( const T& value );
   template< typename T > inline VectorAccessProxy& operator*=( const T& value );
   template< typename T > inline VectorAccessProxy& operator/=( const T& value );
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator Reference() const;
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   Reference get()                       const;
   void      set( ConstReference value ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT&    sv_;  //!< Reference to the accessed sparse vector.
   size_t i_;   //!< Index of the accessed sparse vector element.
   //@}
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
   const Iterator element( sv_.find( i_ ) );
   if( element == sv_.end() )
      sv_.insert( i_, ElementType() );
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
   const Iterator element( sv_.find( i_ ) );
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
inline VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator=( const VectorAccessProxy& vap )
{
   set( vap.get() );
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
inline VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator=( const T& value )
{
   set( value );
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
inline VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator+=( const T& value )
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
inline VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator-=( const T& value )
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
inline VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator*=( const T& value )
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
inline VectorAccessProxy<VT>& VectorAccessProxy<VT>::operator/=( const T& value )
{
   get() /= value;
   return *this;
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
// \return Reference to the accessed sparse vector element.
*/
template< typename VT >  // Type of the sparse vector
inline VectorAccessProxy<VT>::operator Reference() const
{
   return get();
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
// \return Reference to the accessed sparse vector element.
*/
template< typename VT >  // Type of the sparse vector
inline typename VectorAccessProxy<VT>::Reference VectorAccessProxy<VT>::get() const
{
   const Iterator element( sv_.find( i_ ) );
   BLAZE_INTERNAL_ASSERT( element != sv_.end(), "Missing vector element detected" );
   return element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the value of the accessed sparse vector element.
//
// \param value Reference to the new value of the sparse vector element.
// \return void
*/
template< typename VT >  // Type of the sparse vector
inline void VectorAccessProxy<VT>::set( ConstReference value ) const
{
   const Iterator element( sv_.find( i_ ) );
   BLAZE_INTERNAL_ASSERT( element != sv_.end(), "Missing vector element detected" );
   element->value() = value;
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
   typedef typename VectorAccessProxy<VT1>::Reference  LhsReference;
   typedef typename VectorAccessProxy<VT2>::Reference  RhsReference;
   return ( static_cast<LhsReference>( lhs ) == static_cast<RhsReference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( static_cast<Reference>( lhs ) == rhs );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( lhs == static_cast<Reference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT1>::Reference  LhsReference;
   typedef typename VectorAccessProxy<VT2>::Reference  RhsReference;
   return ( static_cast<LhsReference>( lhs ) != static_cast<RhsReference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( static_cast<Reference>( lhs ) != rhs );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( lhs != static_cast<Reference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT1>::Reference  LhsReference;
   typedef typename VectorAccessProxy<VT2>::Reference  RhsReference;
   return ( static_cast<LhsReference>( lhs ) < static_cast<RhsReference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( static_cast<Reference>( lhs ) < rhs );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( lhs < static_cast<Reference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT1>::Reference  LhsReference;
   typedef typename VectorAccessProxy<VT2>::Reference  RhsReference;
   return ( static_cast<LhsReference>( lhs ) > static_cast<RhsReference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( static_cast<Reference>( lhs ) > rhs );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( lhs > static_cast<Reference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT1>::Reference  LhsReference;
   typedef typename VectorAccessProxy<VT2>::Reference  RhsReference;
   return ( static_cast<LhsReference>( lhs ) <= static_cast<RhsReference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( static_cast<Reference>( lhs ) <= rhs );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( lhs <= static_cast<Reference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT1>::Reference  LhsReference;
   typedef typename VectorAccessProxy<VT2>::Reference  RhsReference;
   return ( static_cast<LhsReference>( lhs ) >= static_cast<RhsReference>( rhs ) );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( static_cast<Reference>( lhs ) >= rhs );
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
   typedef typename VectorAccessProxy<VT>::Reference  Reference;
   return ( lhs >= static_cast<Reference>( rhs ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
