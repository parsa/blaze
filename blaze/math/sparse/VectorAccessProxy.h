//=================================================================================================
/*!
//  \file blaze/math/sparse/VectorAccessProxy.h
//  \brief Header file for the VectorAccessProxy class
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

#ifndef _BLAZE_MATH_SPARSE_VECTORACCESSPROXY_H_
#define _BLAZE_MATH_SPARSE_VECTORACCESSPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

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
   // No explicitly declared destructor.
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

} // namespace blaze

#endif
