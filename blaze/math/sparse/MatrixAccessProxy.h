//=================================================================================================
/*!
//  \file blaze/math/sparse/MatrixAccessProxy.h
//  \brief Header file for the MatrixAccessProxy class
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

#ifndef _BLAZE_MATH_SPARSE_MATRIXACCESSPROXY_H_
#define _BLAZE_MATH_SPARSE_MATRIXACCESSPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <ostream>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/proxy/Proxy.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for sparse, \f$ M \times N \f$ matrices.
// \ingroup math
//
// The MatrixAccessProxy provides safe access to the elements of a non-const sparse matrices.\n
// The proxied access to the elements of a sparse matrix is necessary since it may be possible
// that several insertion operations happen in the same statement. The following code illustrates
// this with two examples by means of the CompressedMatrix class:

   \code
   blaze::CompressedMatrix<double> A( 4, 4 );

   // Standard usage of the function call operator to initialize a matrix element.
   // Only a single sparse matrix element is accessed!
   A(0,1) = 1.0;

   // Initialization of a matrix element via another matrix element.
   // Two sparse matrix accesses in one statement!
   A(1,2) = A(0,1);

   // Multiple accesses to elements of the sparse matrix in one statement!
   const double result = A(0,2) + A(1,2) + A(2,2);
   \endcode

// The problem (especially with the last statement) is that several insertion operations might
// take place due to the access via the function call operator. If the function call operator
// would return a direct reference to one of the accessed elements, this reference might be
// invalidated during the evaluation of a subsequent function call operator, which results in
// undefined behavior. This class provides the necessary functionality to guarantee a safe access
// to the sparse matrix elements while preserving the intuitive use of the function call operator.
*/
template< typename MT >  // Type of the sparse matrix
class MatrixAccessProxy : public Proxy< MatrixAccessProxy<MT>, typename MT::ElementType >
{
 private:
   //**Enumerations********************************************************************************
   //! Compile time flag indicating whether the given matrix type is a row-major matrix.
   enum { rmm = IsRowMajorMatrix<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef typename MT::ElementType  RepresentedType;  //!< Type of the represented sparse matrix element.
   typedef RepresentedType&          RawReference;     //!< Raw reference to the represented element.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline MatrixAccessProxy( MT& sm, size_t i, size_t j );
            inline MatrixAccessProxy( const MatrixAccessProxy& map );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~MatrixAccessProxy();
   //@}
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                          inline const MatrixAccessProxy& operator= ( const MatrixAccessProxy& map ) const;
   template< typename T > inline const MatrixAccessProxy& operator= ( const T& value ) const;
   template< typename T > inline const MatrixAccessProxy& operator+=( const T& value ) const;
   template< typename T > inline const MatrixAccessProxy& operator-=( const T& value ) const;
   template< typename T > inline const MatrixAccessProxy& operator*=( const T& value ) const;
   template< typename T > inline const MatrixAccessProxy& operator/=( const T& value ) const;
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
   MT&    sm_;  //!< Reference to the accessed sparse matrix.
   size_t i_;   //!< Row-index of the accessed sparse matrix element.
   size_t j_;   //!< Column-index of the accessed sparse matrix element.
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
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
/*!\brief Initialization constructor for a MatrixAccessProxy.
//
// \param sm Reference to the accessed sparse matrix.
// \param i The row-index of the accessed sparse matrix element.
// \param j The column-index of the accessed sparse matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::MatrixAccessProxy( MT& sm, size_t i, size_t j )
   : sm_( sm )  // Reference to the accessed sparse matrix
   , i_ ( i  )  // Row-index of the accessed sparse matrix element
   , j_ ( j  )  // Column-index of the accessed sparse matrix element
{
   const typename MT::Iterator element( sm_.find( i_, j_ ) );
   const size_t index( rmm ? i_ : j_ );
   if( element == sm_.end(index) )
      sm_.insert( i_, j_, RepresentedType() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for MatrixAccessProxy.
//
// \param map Sparse matrix access proxy to be copied.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::MatrixAccessProxy( const MatrixAccessProxy& map )
   : sm_( map.sm_ )  // Reference to the accessed sparse matrix
   , i_ ( map.i_  )  // Row-index of the accessed sparse matrix element
   , j_ ( map.j_  )  // Column-index of the accessed sparse matrix element
{
   BLAZE_INTERNAL_ASSERT( sm_.find(i_,j_) != sm_.end( rmm ? i_ : j_ ), "Missing matrix element detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for MatrixAccessProxy.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::~MatrixAccessProxy()
{
   const typename MT::Iterator element( sm_.find( i_, j_ ) );
   const size_t index( rmm ? i_ : j_ );
   if( element != sm_.end( index ) && isDefault( element->value() ) )
      sm_.erase( index, element );
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for MatrixAccessProxy.
//
// \param map Sparse matrix access proxy to be copied.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
inline const MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator=( const MatrixAccessProxy& map ) const
{
   get() = map.get();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the accessed sparse matrix element.
//
// \param value The new value of the sparse matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline const MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator=( const T& value ) const
{
   get() = value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value to be added to the sparse matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline const MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator+=( const T& value ) const
{
   get() += value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value to be subtracted from the sparse matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline const MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator-=( const T& value ) const
{
   get() -= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline const MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator*=( const T& value ) const
{
   get() *= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline const MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator/=( const T& value ) const
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
/*!\brief Returning the value of the accessed sparse matrix element.
//
// \return Direct/raw reference to the accessed sparse matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::RawReference MatrixAccessProxy<MT>::get() const
{
   const typename MT::Iterator element( sm_.find( i_, j_ ) );
   BLAZE_INTERNAL_ASSERT( element != sm_.end( rmm ? i_ : j_ ), "Missing matrix element detected" );
   return element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the proxy represents a restricted sparse matrix element..
//
// \return \a true in case access to the sparse matrix element is restricted, \a false if not.
*/
template< typename MT >  // Type of the sparse matrix
inline bool MatrixAccessProxy<MT>::isRestricted() const
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
/*!\brief Conversion to the accessed sparse matrix element.
//
// \return Direct/raw reference to the accessed sparse matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::operator RawReference() const
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
/*!\name MatrixAccessProxy operators */
//@{
template< typename MT1, typename MT2 >
inline bool operator==( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator==( const MatrixAccessProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator==( const T& lhs, const MatrixAccessProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator!=( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator!=( const MatrixAccessProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator!=( const T& lhs, const MatrixAccessProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<( const MatrixAccessProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<( const T& lhs, const MatrixAccessProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>( const MatrixAccessProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>( const T& lhs, const MatrixAccessProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<=( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<=( const MatrixAccessProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<=( const T& lhs, const MatrixAccessProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>=( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>=( const MatrixAccessProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>=( const T& lhs, const MatrixAccessProxy<MT>& rhs );

template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const MatrixAccessProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two MatrixAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename MT1, typename MT2 >
inline bool operator==( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs )
{
   return ( lhs.get() == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a MatrixAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename MT, typename T >
inline bool operator==( const MatrixAccessProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a MatrixAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename MT >
inline bool operator==( const T& lhs, const MatrixAccessProxy<MT>& rhs )
{
   return ( lhs == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two MatrixAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename MT1, typename MT2 >
inline bool operator!=( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs )
{
   return ( lhs.get() != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a MatrixAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename MT, typename T >
inline bool operator!=( const MatrixAccessProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a MatrixAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename MT >
inline bool operator!=( const T& lhs, const MatrixAccessProxy<MT>& rhs )
{
   return ( lhs != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two MatrixAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs )
{
   return ( lhs.get() < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a MatrixAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<( const MatrixAccessProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a MatrixAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<( const T& lhs, const MatrixAccessProxy<MT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two MatrixAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs )
{
   return ( lhs.get() > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a MatrixAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>( const MatrixAccessProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a MatrixAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>( const T& lhs, const MatrixAccessProxy<MT>& rhs )
{
   return ( lhs > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two MatrixAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<=( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs )
{
   return ( lhs.get() <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a MatrixAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<=( const MatrixAccessProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a MatrixAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<=( const T& lhs, const MatrixAccessProxy<MT>& rhs )
{
   return ( lhs <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two MatrixAccessProxy objects.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>=( const MatrixAccessProxy<MT1>& lhs, const MatrixAccessProxy<MT2>& rhs )
{
   return ( lhs.get() >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a MatrixAccessProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side MatrixAccessProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>=( const MatrixAccessProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a MatrixAccessProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side MatrixAccessProxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>=( const T& lhs, const MatrixAccessProxy<MT>& rhs )
{
   return ( lhs >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the MatrixAccessProxy class template.
// \ingroup math
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant proxy object.
// \return Reference to the output stream.
*/
template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const MatrixAccessProxy<MT>& proxy )
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
/*!\name MatrixAccessProxy global functions */
//@{
template< typename MT >
inline void reset( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline void clear( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline bool isDefault( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline void swap( const MatrixAccessProxy<MT>& a, const MatrixAccessProxy<MT>& b ) /* throw() */;

template< typename MT, typename T >
inline void swap( const MatrixAccessProxy<MT>& a, T& b ) /* throw() */;

template< typename T, typename MT >
inline void swap( T& a, const MatrixAccessProxy<MT>& v ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the represented element to the default initial values.
// \ingroup math
//
// \param proxy The given access proxy.
// \return void
//
// This function resets the element represented by the access proxy to its default initial value.
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// reset() function, this function resets all elements of the vector/matrix to the default initial
// values.
*/
template< typename MT >
inline void reset( const MatrixAccessProxy<MT>& proxy )
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
//
// This function clears the element represented by the access proxy to its default initial state.
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// clear() function, this function clears the vector/matrix to its default initial state.
*/
template< typename MT >
inline void clear( const MatrixAccessProxy<MT>& proxy )
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
template< typename MT >
inline bool isDefault( const MatrixAccessProxy<MT>& proxy )
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
template< typename MT >
inline void swap( const MatrixAccessProxy<MT>& a, const MatrixAccessProxy<MT>& b ) /* throw() */
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
template< typename MT, typename T >
inline void swap( const MatrixAccessProxy<MT>& a, T& b ) /* throw() */
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
template< typename T, typename MT >
inline void swap( T& a, const MatrixAccessProxy<MT>& b ) /* throw() */
{
   using std::swap;

   swap( a, b.get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
