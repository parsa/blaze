//=================================================================================================
/*!
//  \file blaze/math/adaptors/symmetricmatrix/NonNumericProxy.h
//  \brief Header file for the NonNumericProxy class
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_NONNUMERICPROXY_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_NONNUMERICPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <ostream>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/proxy/Proxy.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for symmetric, square matrices with non-numeric element types.
// \ingroup symmetric_matrix
//
// The NonNumericProxy provides controlled access to the elements of a non-const symmetric matrix
// with non-numeric element type (e.g. vectors or matrices). It guarantees that a modification of
// element \f$ a_{ij} \f$ of the accessed matrix is also applied to element \f$ a_{ji} \f$. The
// following example illustrates this by means of a \f$ 3 \times 3 \f$ sparse symmetric matrix
// with StaticVector elements:

   \code
   using blaze::CompressedMatrix;
   using blaze::StaticVector;
   using blaze::SymmetricMatrix;

   typedef StaticVector<int,3UL>  Vector;

   // Creating a 3x3 symmetric sparses matrix
   SymmetricMatrix< CompressedMatrix< Vector > > A( 3UL );

   A(0,2) = Vector( -2,  1 );  //        ( (  0 0 ) ( 0  0 ) ( -2  1 ) )
   A(1,1) = Vector(  3,  4 );  // => A = ( (  0 0 ) ( 3  4 ) (  5 -1 ) )
   A(1,2) = Vector(  5, -1 );  //        ( ( -2 1 ) ( 5 -1 ) (  0  0 ) )
   \endcode
*/
template< typename MT >  // Type of the adapted matrix
class NonNumericProxy : public Proxy< NonNumericProxy<MT>, typename MT::ElementType::ValueType >
{
 private:
   //**Enumerations********************************************************************************
   //! Compile time flag indicating whether the given matrix type is a row-major matrix.
   enum { rmm = IsRowMajorMatrix<MT>::value };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename MT::ElementType  ET;  //!< Element type of the adapted matrix.
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef typename ET::ValueType  RepresentedType;  //!< Type of the represented matrix element.
   typedef typename ET::Reference  RawReference;     //!< Raw reference to the represented element.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline NonNumericProxy( MT& sm, size_t i, size_t j );
            inline NonNumericProxy( const NonNumericProxy& nnp );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~NonNumericProxy();
   //@}
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                          inline NonNumericProxy& operator= ( const NonNumericProxy& nnp );
   template< typename T > inline NonNumericProxy& operator= ( const T& value );
   template< typename T > inline NonNumericProxy& operator+=( const T& value );
   template< typename T > inline NonNumericProxy& operator-=( const T& value );
   template< typename T > inline NonNumericProxy& operator*=( const T& value );
   template< typename T > inline NonNumericProxy& operator/=( const T& value );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline RawReference get() const;
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
   MT&    matrix_;  //!< Reference to the adapted matrix.
   size_t i_;       //!< Row-index of the accessed matrix element.
   size_t j_;       //!< Column-index of the accessed matrix element.
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST                ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE             ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_NUMERIC_TYPE         ( RepresentedType );
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
/*!\brief Initialization constructor for a NonNumericProxy.
//
// \param matrix Reference to the adapted matrix.
// \param i The row-index of the accessed matrix element.
// \param j The column-index of the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline NonNumericProxy<MT>::NonNumericProxy( MT& matrix, size_t i, size_t j )
   : matrix_( matrix )  // Reference to the adapted matrix
   , i_     ( i )       // Row-index of the accessed matrix element
   , j_     ( j )       // Column-index of the accessed matrix element
{
   const typename MT::Iterator pos( matrix_.find( i_, j_ ) );
   const size_t index( rmm ? i_ : j_ );

   if( pos == matrix_.end(index) )
   {
      const typename MT::ElementType element( ( RepresentedType() ) );
      matrix_.insert( i_, j_, element );
      if( i_ != j_ )
         matrix_.insert( j_, i_, element );
   }

   BLAZE_INTERNAL_ASSERT( matrix_.find(i_,j_)->value() == matrix_.find(j_,i_)->value(), "Unbalance detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for NonNumericProxy.
//
// \param nnp Non-numeric access proxy to be copied.
*/
template< typename MT >  // Type of the adapted matrix
inline NonNumericProxy<MT>::NonNumericProxy( const NonNumericProxy& nnp )
   : matrix_( nnp.matrix_ )  // Reference to the adapted matrix
   , i_     ( nnp.i_ )       // Row-index of the accessed matrix element
   , j_     ( nnp.j_ )       // Column-index of the accessed matrix element
{
   BLAZE_INTERNAL_ASSERT( matrix_.find(i_,j_) != matrix_.end( rmm ? i_ : j_ ), "Missing matrix element detected" );
   BLAZE_INTERNAL_ASSERT( matrix_.find(j_,i_) != matrix_.end( rmm ? j_ : i_ ), "Missing matrix element detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for NonNumericProxy.
*/
template< typename MT >  // Type of the adapted matrix
inline NonNumericProxy<MT>::~NonNumericProxy()
{
   const typename MT::Iterator pos( matrix_.find( i_, j_ ) );
   const size_t index( rmm ? i_ : j_ );

   if( pos != matrix_.end( index ) && isDefault( *pos->value() ) )
   {
      matrix_.erase( index, pos );
      if( i_ != j_ )
         matrix_.erase( ( rmm ? j_ : i_ ), matrix_.find( j_, i_ ) );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for NonNumericProxy.
//
// \param nnp Non-numeric access proxy to be copied.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the adapted matrix
inline NonNumericProxy<MT>& NonNumericProxy<MT>::operator=( const NonNumericProxy& nnp )
{
   get() = nnp.get();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the represented matrix element.
//
// \param value The new value of the matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NonNumericProxy<MT>& NonNumericProxy<MT>::operator=( const T& value )
{
   get() = value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the represented matrix element.
//
// \param value The right-hand side value to be added to the matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NonNumericProxy<MT>& NonNumericProxy<MT>::operator+=( const T& value )
{
   get() += value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the represented matrix element.
//
// \param value The right-hand side value to be subtracted from the matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NonNumericProxy<MT>& NonNumericProxy<MT>::operator-=( const T& value )
{
   get() -= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the represented matrix element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NonNumericProxy<MT>& NonNumericProxy<MT>::operator*=( const T& value )
{
   get() *= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the represented matrix element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NonNumericProxy<MT>& NonNumericProxy<MT>::operator/=( const T& value )
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
/*!\brief Returning a reference to the accessed matrix element.
//
// \return Direct/raw reference to the accessed matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline typename NonNumericProxy<MT>::RawReference NonNumericProxy<MT>::get() const
{
   const typename MT::Iterator pos( matrix_.find( i_, j_ ) );
   BLAZE_INTERNAL_ASSERT( pos != matrix_.end( rmm ? i_ : j_ ), "Missing matrix element detected" );
   return *pos->value();
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion to the represented matrix element.
//
// \return Direct/raw reference to the represented matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline NonNumericProxy<MT>::operator RawReference() const
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
/*!\name NonNumericProxy operators */
//@{
template< typename MT1, typename MT2 >
inline bool operator==( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator==( const NonNumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator==( const T& lhs, const NonNumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator!=( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator!=( const NonNumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator!=( const T& lhs, const NonNumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<( const NonNumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<( const T& lhs, const NonNumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>( const NonNumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>( const T& lhs, const NonNumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<=( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<=( const NonNumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<=( const T& lhs, const NonNumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>=( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>=( const NonNumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>=( const T& lhs, const NonNumericProxy<MT>& rhs );

template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const NonNumericProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two NonNumericProxy objects.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename MT1, typename MT2 >
inline bool operator==( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs )
{
   return ( lhs.get() == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a NonNumericProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename MT, typename T >
inline bool operator==( const NonNumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a NonNumericProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename MT >
inline bool operator==( const T& lhs, const NonNumericProxy<MT>& rhs )
{
   return ( lhs == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two NonNumericProxy objects.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename MT1, typename MT2 >
inline bool operator!=( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs )
{
   return ( lhs.get() != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a NonNumericProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename MT, typename T >
inline bool operator!=( const NonNumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a NonNumericProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename MT >
inline bool operator!=( const T& lhs, const NonNumericProxy<MT>& rhs )
{
   return ( lhs != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two NonNumericProxy objects.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs )
{
   return ( lhs.get() < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a NonNumericProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<( const NonNumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a NonNumericProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<( const T& lhs, const NonNumericProxy<MT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two NonNumericProxy objects.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs )
{
   return ( lhs.get() > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a NonNumericProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>( const NonNumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a NonNumericProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>( const T& lhs, const NonNumericProxy<MT>& rhs )
{
   return ( lhs > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two NonNumericProxy objects.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<=( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs )
{
   return ( lhs.get() <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a NonNumericProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<=( const NonNumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a NonNumericProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<=( const T& lhs, const NonNumericProxy<MT>& rhs )
{
   return ( lhs <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two NonNumericProxy objects.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>=( const NonNumericProxy<MT1>& lhs, const NonNumericProxy<MT2>& rhs )
{
   return ( lhs.get() >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a NonNumericProxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side NonNumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>=( const NonNumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a NonNumericProxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NonNumericProxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>=( const T& lhs, const NonNumericProxy<MT>& rhs )
{
   return ( lhs >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the NonNumericProxy class template.
// \ingroup math
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant proxy object.
// \return Reference to the output stream.
*/
template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const NonNumericProxy<MT>& proxy )
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
/*!\name NonNumericProxy global functions */
//@{
template< typename MT >
BLAZE_ALWAYS_INLINE void reset( const NonNumericProxy<MT>& proxy );

template< typename MT >
BLAZE_ALWAYS_INLINE void clear( const NonNumericProxy<MT>& proxy );

template< typename MT >
BLAZE_ALWAYS_INLINE bool isDefault( const NonNumericProxy<MT>& proxy );
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
BLAZE_ALWAYS_INLINE void reset( const NonNumericProxy<MT>& proxy )
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
BLAZE_ALWAYS_INLINE void clear( const NonNumericProxy<MT>& proxy )
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
BLAZE_ALWAYS_INLINE bool isDefault( const NonNumericProxy<MT>& proxy )
{
   using blaze::isDefault;

   return isDefault( proxy.get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
