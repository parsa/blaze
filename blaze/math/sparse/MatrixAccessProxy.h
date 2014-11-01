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

#include <ostream>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/GetMemberType.h>


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
class MatrixAccessProxy
{
 private:
   //**Enumerations********************************************************************************
   //! Compile time flag indicating whether the given matrix type is a row-major matrix.
   enum { rmm = IsRowMajorMatrix<MT>::value };
   //**********************************************************************************************

   //**Type trait generation***********************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetValueType    , value_type   , INVALID_TYPE );
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetElementType  , ElementType  , INVALID_TYPE );
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetReference    , Reference    , INVALID_TYPE );
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetPointer      , Pointer      , INVALID_TYPE );
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetIterator     , Iterator     , INVALID_TYPE );
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetConstIterator, ConstIterator, INVALID_TYPE );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef typename MT::ElementType  RepresentedType;  //!< Type of the represented sparse matrix element.
   typedef RepresentedType&          RawReference;     //!< Raw reference to the represented element.

   //! Value type of the represented complex element.
   typedef typename GetValueType<RepresentedType>::Type  ValueType;

   //! Element type of the represented sparse vector element.
   typedef typename GetElementType<RepresentedType>::Type  ElementType;

   //! Reference type of the represented sparse matrix element.
   typedef typename GetReference<RepresentedType>::Type  Reference;

   //! Pointer type of the represented sparse matrix element.
   typedef typename GetPointer<RepresentedType>::Type  Pointer;

   //! Iterator type of the represented sparse matrix element.
   typedef typename GetIterator<RepresentedType>::Type  Iterator;

   //! ConstIterator type of the represented sparse matrix element.
   typedef typename GetConstIterator<RepresentedType>::Type  ConstIterator;
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
                          inline MatrixAccessProxy& operator= ( const MatrixAccessProxy& map );
   template< typename T > inline MatrixAccessProxy& operator= ( const T& value );
   template< typename T > inline MatrixAccessProxy& operator+=( const T& value );
   template< typename T > inline MatrixAccessProxy& operator-=( const T& value );
   template< typename T > inline MatrixAccessProxy& operator*=( const T& value );
   template< typename T > inline MatrixAccessProxy& operator/=( const T& value );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void reset() const;
   inline void clear() const;
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator RawReference() const;
   //@}
   //**********************************************************************************************

   //**Complex data access functions***************************************************************
   /*!\name Complex data access functions */
   //@{
   inline ValueType real() const;
   inline void      real( ValueType value ) const;
   inline ValueType imag() const;
   inline void      imag( ValueType value ) const;
   //@}
   //**********************************************************************************************

   //**Vector/matrix data access functions*********************************************************
   /*!\name Vector/matrix data access functions */
   //@{
   inline Reference operator[]( size_t index ) const;
   inline Reference operator()( size_t i, size_t j ) const;

   inline Pointer       data  () const;
   inline Pointer       data  ( size_t i ) const;
   inline Iterator      begin () const;
   inline Iterator      begin ( size_t i ) const;
   inline ConstIterator cbegin() const;
   inline ConstIterator cbegin( size_t i ) const;
   inline Iterator      end   () const;
   inline Iterator      end   ( size_t i ) const;
   inline ConstIterator cend  () const;
   inline ConstIterator cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Vector/matrix utility functions*************************************************************
   /*!\name Vector/matrix utility functions */
   //@{
   inline size_t   size() const;
   inline size_t   rows() const;
   inline size_t   columns() const;
   inline size_t   spacing() const;
   inline size_t   capacity() const;
   inline size_t   capacity( size_t i ) const;
   inline size_t   nonZeros() const;
   inline size_t   nonZeros( size_t i ) const;
   inline void     reset( size_t i ) const;
   inline Iterator set( size_t index, const ElementType& value ) const;
   inline Iterator set( size_t i, size_t j, const ElementType& value ) const;
   inline Iterator insert( size_t index, const ElementType& value ) const;
   inline Iterator insert( size_t i, size_t j, const ElementType& value ) const;
   inline void     append( size_t index, const ElementType& value, bool check=false ) const;
   inline void     append( size_t i, size_t j, const ElementType& value, bool check=false ) const;
   inline void     finalize( size_t i ) const;
   inline void     erase( size_t index ) const;
   inline Iterator erase( Iterator pos ) const;
   inline Iterator erase( Iterator first, Iterator last ) const;
   inline void     erase( size_t i, size_t j ) const;
   inline Iterator erase( size_t i, Iterator pos ) const;
   inline Iterator erase( size_t i, Iterator first, Iterator last ) const;
   inline void     resize( size_t n, bool preserve=true ) const;
   inline void     resize ( size_t m, size_t n, bool preserve=true ) const;
   inline void     extend( size_t n, bool preserve=true ) const;
   inline void     extend ( size_t m, size_t n, bool preserve=true ) const;
   inline void     reserve( size_t n ) const;
   inline void     reserve( size_t i, size_t n ) const;
   inline void     trim() const;
   inline void     trim( size_t i ) const;
   inline void     transpose() const;

   template< typename Other > inline void scale( const Other& scalar ) const;
   //@}
   //**********************************************************************************************

   //**Vector/matrix lookup functions**************************************************************
   /*!\name Vector/matrix lookup functions */
   //@{
   inline Iterator find      ( size_t index ) const;
   inline Iterator find      ( size_t i, size_t j ) const;
   inline Iterator lowerBound( size_t index ) const;
   inline Iterator lowerBound( size_t i, size_t j ) const;
   inline Iterator upperBound( size_t index ) const;
   inline Iterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   RawReference get() const;
   //@}
   //**********************************************************************************************

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
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator=( const MatrixAccessProxy& map )
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
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator=( const T& value )
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
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator+=( const T& value )
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
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator-=( const T& value )
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
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator*=( const T& value )
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
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator/=( const T& value )
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
/*!\brief Reset to the default initial values.
//
// \return void
//
// This function resets the element represented by the access proxy to its default initial value.
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// reset() function, this function resets all elements of the vector/matrix to the default initial
// values.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::reset() const
{
   using blaze::reset;

   reset( get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
//
// \return void
//
// This function clears the element represented by the access proxy to its default initial state.
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// clear() function, this function clears the vector/matrix to its default initial state.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::clear() const
{
   using blaze::clear;

   clear( get() );
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
//  COMPLEX DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the real part of the represented complex number.
//
// \return The current real part of the represented complex number.
//
// In case the access proxy represents a complex number, this function returns the current value
// of its real part.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::ValueType MatrixAccessProxy<MT>::real() const
{
   return get().real();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the real part of the represented complex number.
//
// \param value The new value for the real part.
// \return void
//
// In case the access proxy represents a complex number, this function sets a new value to its
// real part.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::real( ValueType value ) const
{
   get().real( value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the imaginary part of the represented complex number.
//
// \return The current imaginary part of the represented complex number.
//
// In case the access proxy represents a complex number, this function returns the current value
// of its imaginary part.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::ValueType MatrixAccessProxy<MT>::imag() const
{
   return get().imag();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the imaginary part of the represented complex number.
//
// \param value The new value for the imaginary part.
// \return void
//
// In case the access proxy represents a complex number, this function sets a new value to its
// imaginary part.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::imag( ValueType value ) const
{
   get().imag( value );
}
//*************************************************************************************************




//=================================================================================================
//
//  VECTOR/MATRIX UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// In case the access proxy represents a vector-like data structure that provides a subscript
// operator, this function provides access to the vector elements via this subscript operator.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Reference
   MatrixAccessProxy<MT>::operator[]( size_t index ) const
{
   return get()[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Function call operator for the direct access to matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// In case the access proxy represents a matrix-like data structure that provides a function call
// operator, this function provides access to the matrix elements via this function call operator.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Reference
   MatrixAccessProxy<MT>::operator()( size_t i, size_t j ) const
{
   return get()(i,j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to vector/matrix elements.
//
// \return Pointer to the internal element storage.
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// data() function, this function returns a pointer to the internal storage of the vector/matrix.
// Note that in case the proxy represents a matrix-like data structure you can NOT assume that
// all matrix elements lie adjacent to each other! The matrix may use techniques such as padding
// to improve the alignment of the data. Whereas the number of elements within a row/column are
// given by the \c rows() and \c columns() member functions, respectively, the total number of
// elements including padding is given by the \c spacing() member function.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Pointer MatrixAccessProxy<MT>::data() const
{
   return get().data();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to matrix elements of row/column \a i.
//
// \return Pointer to the internal element storage.
//
// In case the access proxy represents a matrix-like data structure that provides a data()
// function, this function returns a pointer to the internal storage for the elements in
// row/column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Pointer MatrixAccessProxy<MT>::data( size_t i ) const
{
   return get().data(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the represented vector.
//
// \return Iterator to the first element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a begin()
// function, this function returns an iterator to the first element of the vector.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator MatrixAccessProxy<MT>::begin() const
{
   return get().begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i of the represented matrix.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a begin()
// function, this function returns a row/column iterator to the first element of row/column
// \a i. In case the storage order of the matrix is set to \a rowMajor the function returns
// an iterator to the first element of row \a i, in case the storage flag is set to
// \a columnMajor the function returns an iterator to the first element of column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator MatrixAccessProxy<MT>::begin( size_t i ) const
{
   return get().begin(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the represented vector.
//
// \return Iterator to the first element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a cbegin()
// function, this function returns an iterator to the first element of the vector.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::ConstIterator MatrixAccessProxy<MT>::cbegin() const
{
   return get().cbegin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i of the represented matrix.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a cbegin()
// function, this function returns a row/column iterator to the first element of row/column
// \a i. In case the storage order of the matrix is set to \a rowMajor the function returns
// an iterator to the first element of row \a i, in case the storage flag is set to
// \a columnMajor the function returns an iterator to the first element of column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::ConstIterator MatrixAccessProxy<MT>::cbegin( size_t i ) const
{
   return get().cbegin(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the represented vector.
//
// \return Iterator just past the last element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides an end()
// function, this function returns an iterator just past the last element of the vector.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator MatrixAccessProxy<MT>::end() const
{
   return get().end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i of the represented matrix.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides an end()
// function, this function returns an row/column iterator just past the last element of row/column
// \a i. In case the storage order of the matrix is set to \a rowMajor the function returns
// an iterator just past the last element of row \a i, in case the storage flag is set to
// \a columnMajor the function returns an iterator just past the last element of column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator MatrixAccessProxy<MT>::end( size_t i ) const
{
   return get().end(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the represented vector.
//
// \return Iterator just past the last element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a cend()
// function, this function returns an iterator just past the last element of the vector.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::ConstIterator MatrixAccessProxy<MT>::cend() const
{
   return get().cend();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i of the represented matrix.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a cend()
// function, this function returns an row/column iterator just past the last element of row/column
// \a i. In case the storage order of the matrix is set to \a rowMajor the function returns
// an iterator just past the last element of row \a i, in case the storage flag is set to
// \a columnMajor the function returns an iterator just past the last element of column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::ConstIterator MatrixAccessProxy<MT>::cend( size_t i ) const
{
   return get().cend(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size/dimension of the represented vector.
//
// \return The size of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a size()
// function, this function returns the current size/dimension of the vector.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::size() const
{
   return get().size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of rows of the represented matrix.
//
// \return The number of rows of the matrix.
//
// In case the access proxy represents a matrix-like data structure that provides a rows()
// function, this function returns the current number of rows of the matrix.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::rows() const
{
   return get().rows();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the represented matrix.
//
// \return The number of columns of the matrix.
//
// In case the access proxy represents a matrix-like data structure that provides a columns()
// function, this function returns the current number of columns of the matrix.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::columns() const
{
   return get().columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows/columns of the represented matrix.
//
// \return The spacing between the beginning of two rows/columns.
//
// In case the access proxy represents a matrix-like data structure that provides a spacing()
// function, this function returns the spacing between the beginning of two rows/columns, i.e.
// the total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::spacing() const
{
   return get().spacing();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the represented vector/matrix.
//
// \return The capacity of the vector/matrix.
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// capacity() function, this function returns the current maximum capacity of the vector/matrix.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::capacity() const
{
   return get().capacity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row/column of the represented matrix.
//
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a capacity()
// function, this function returns the current capacity of the specified row/column. In case
// the storage order is set to \a rowMajor the function returns the capacity of row \a i, in
// case the storage flag is set to \a columnMajor the function returns the capacity of column
// \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::capacity( size_t i ) const
{
   return get().capacity(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the represented vector/matrix.
//
// \return The number of non-zero elements in the vector/matrix.
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// nonZeros() function, this function returns the current number of non-zero elements in the
// vector/matrix.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::nonZeros() const
{
   return get().nonZeros();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column of the represented matrix.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a nonZeros()
// function, this function returns the current number of non-zero elements in the specified
// row/column. In case the storage order is set to \a rowMajor the function returns the number
// of non-zero elements in row \a i, in case the storage flag is set to \a columnMajor the
// function returns the number of non-zero elements in column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t MatrixAccessProxy<MT>::nonZeros( size_t i ) const
{
   return get().nonZeros(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column to the default initial values.
//
// \param i The index of the row/column.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a reset()
// function, this function resets all elements of the matrix to the default initial values.
// In case the storage order is set to \a rowMajor the function resets the values in row \a i,
// in case the storage order is set to \a columnMajor the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::reset( size_t i ) const
{
   using blaze::reset;

   reset( get(), i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an element of the represented sparse vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
// \exception std::invalid_argument Invalid compressed vector access index.
//
// In case the access proxy represents a vector-like data structure that provides a set() function,
// this function sets the value of an element of the sparse vector. In case the sparse vector
// already contains an element with index \a index its value is modified, else a new element with
// the given \a value is inserted.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::set( size_t index, const ElementType& value ) const
{
   return get().set( index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an element of the represented sparse matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
// \exception std::invalid_argument Invalid sparse matrix access index.
//
// In case the access proxy represents a matrix-like data structure that provides a set() function,
// this function sets the value of an element of the sparse matrix. In case the sparse matrix
// already contains an element with row index \a i and column index \a j its value is modified,
// else a new element with the given \a value is inserted.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::set( size_t i, size_t j, const ElementType& value ) const
{
   return get().set( i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the represented sparse vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid compressed vector access index.
//
// In case the access proxy represents a vector-like data structure that provides an insert()
// function, this function inserts a new element into the sparse vector. However, duplicate
// elements are not allowed. In case the sparse vector already contains an element with index
// \a index, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::insert( size_t index, const ElementType& value ) const
{
   return get().insert( index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the represented sparse matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid sparse matrix access index.
//
// In case the access proxy represents a vector-like data structure that provides an insert()
// function, this function inserts a new element into the sparse matrix. However, duplicate
// elements are not allowed. In case the sparse matrix already contains an element with row
// index \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::insert( size_t i, size_t j, const ElementType& value ) const
{
   return get().insert( i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Appending an element to the represented sparse vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// In case the access proxy represents a vector-like data structure that provides an append()
// function, this function provides a very efficient way to fill a sparse vector with elements.
// It appends a new element to the end of the sparse vector without any memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the compressed vector
//  - the current number of non-zero elements must be smaller than the capacity of the vector
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::append( size_t index, const ElementType& value, bool check ) const
{
   get().append( index, value, check );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Appending an element to the specified row/column of the sparse matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse matrix with elements. It appends
// a new element to the end of the specified row/column without any additional memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row/column of the sparse matrix
//  - the current number of non-zero elements in the matrix must be smaller than the capacity
//    of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::append( size_t i, size_t j, const ElementType& value, bool check ) const
{
   get().append( i, j, value, check );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finalizing the element insertion of a row/column.
//
// \param i The index of the row/column to be finalized \f$[0..M-1]\f$.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a finalize()
// function, this function provides the functionality to finalize the element insertion of a
// row/column. This function is part of the low-level interface to efficiently fill a sparse
// matrix with elements. After completion of row/column \a i via the append() function, this
// function can be called to finalize row/column \a i and prepare the next row/column for
// insertion process via append().
//
// \b Note: Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::finalize( size_t i ) const
{
   get().finalize( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the compressed vector.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// In case the access proxy represents a vector-like data structure that provides an erase()
// function, this function erases an element from the sparse vector.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::erase( size_t index ) const
{
   get().erase( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the compressed vector.
//
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// In case the access proxy represents a vector-like data structure that provides an erase()
// function, this function erases an element from the sparse vector.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator MatrixAccessProxy<MT>::erase( Iterator pos ) const
{
   return get().erase( pos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the compressed vector.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// In case the access proxy represents a vector-like data structure that provides an erase()
// function, this function erases a range of elements from the sparse vector.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::erase( Iterator first, Iterator last ) const
{
   return get().erase( first, last );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse matrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides an erase()
// function, this function erases an element from the sparse matrix.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::erase( size_t i, size_t j ) const
{
   get().erase( i, j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse matrix.
//
// \param i The row/column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// In case the access proxy represents a matrix-like data structure that provides an erase()
// function, this function erases an element from the sparse matrix. In case the storage order
// is set to \a rowMajor the function erases an element from row \a i, in case the storage flag
// is set to \a columnMajor the function erases an element from column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::erase( size_t i, Iterator pos ) const
{
   return get().erase( i, pos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the sparse matrix.
//
// \param i The row/column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// In case the access proxy represents a matrix-like data structure that provides an erase()
// function, this function erases a range of element from the sparse matrix. In case the storage
// order is set to \a rowMajor the function erases a range of elements from row \a i, in case the
// storage flag is set to \a columnMajor the function erases a range of elements from column \a i.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::erase( size_t i, Iterator first, Iterator last ) const
{
   return get().erase( i, first, last );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the represented vector.
//
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// In case the access proxy represents a vector-like data structure that provides a resize()
// function, this function changes the size of the vector. Depending on the type of the vector,
// during this operation new dynamic memory may be allocated in case the capacity of the vector
// is too small. Note that this function may invalidate all existing views (subvectors, ...) on
// the vector if it is used to shrink the vector. Additionally, the resize() operation potentially
// changes all vector elements. In order to preserve the old vector values, the \a preserve flag
// can be set to \a true. However, note that depending on the type of the vector new vector
// elements may not initialized!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::resize( size_t n, bool preserve ) const
{
   get().resize( n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the represented matrix.
//
// \param m The new number of rows of the matrix.
// \param n The new number of columns of the matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a resize()
// function, this function resizes the matrix using the given size to \f$ m \times n \f$.
// Depending on the type of the matrix, during this operation new dynamic memory may be allocated
// in case the capacity of the matrix is too small. Note that this function may invalidate all
// existing views (submatrices, rows, columns, ...) on the matrix if it is used to shrink the
// matrix. Additionally, the resize operation potentially changes all matrix elements. In order
// to preserve the old matrix values, the \a preserve flag can be set to \a true. However, note
// that depending on the type of the matrix new matrix elements may not initialized!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::resize( size_t m, size_t n, bool preserve ) const
{
   get().resize( m, n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the represented vector.
//
// \param n Number of additional vector elements.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// In case the access proxy represents a vector-like data structure that provides a extend()
// function, this function extends the size of the vector. Depending on the type of the vector,
// during this operation new dynamic memory may be allocated in case the capacity of the vector
// is too small. Therefore this function potentially changes all vector elements. In order to
// preserve the old vector values, the \a preserve flag can be set to \a true. However, note
// that depending on the type vector new vector elements may not initialized!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::extend( size_t n, bool preserve ) const
{
   get().extend( n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the represented matrix.
//
// \param m Number of additional rows.
// \param n Number of additional columns.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides an extend()
// function, this function increases the matrix size by \a m rows and \a n columns. Depending
// on the type of the matrix, during this operation new dynamic memory may be allocated in case
// the capacity of the matrix is too small. Therefore this function potentially changes all
// matrix elements. In order to preserve the old matrix values, the \a preserve flag can be
// set to \a true. However, note that depending on the type of the matrix new matrix elements
// may not initialized!
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::extend( size_t m, size_t n, bool preserve ) const
{
   get().extend( m, n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the represented vector/matrix.
//
// \param n The new minimum capacity of the vector/matrix.
// \return void
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// reserve() function, this function increases the capacity of the vector/matrix to at least
// \a n elements. The current values of the vector/matrix elements are preserved.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::reserve( size_t n ) const
{
   get().reserve( n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of a specific row/column of the sparse matrix.
//
// \param i The row/column index of the new element \f$[0..M-1]\f$ or \f$[0..N-1]\f$.
// \param n The new minimum capacity of the specified row/column.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a reserve()
// function, this function increases the capacity of row/column \a i of the sparse matrix to
// at least \a n elements. The current values of the sparse matrix and all other individual
// row/column capacities are preserved. In case the storage order is set to \a rowMajor, the
// function reserves capacity for row \a i and the index has to be in the range \f$[0..M-1]\f$.
// In case the storage order is set to \a columnMajor, the function reserves capacity for column
// \a i and the index has to be in the range \f$[0..N-1]\f$.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::reserve( size_t i, size_t n ) const
{
   get().reserve( i, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all excessive capacity from all rows/columns.
//
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a trim()
// function, this function can be used to reverse the effect of all row/column-specific reserve()
// calls. The function removes all excessive capacity from all rows (in case of a rowMajor
// matrix) or columns (in case of a columnMajor matrix). Note that this function does not
// remove the overall capacity but only reduces the capacity per row/column.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::trim() const
{
   get().trim();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all excessive capacity of a specific row/column of the sparse matrix.
//
// \param i The index of the row/column to be trimmed (\f$[0..M-1]\f$ or \f$[0..N-1]\f$).
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a trim()
// function, this function can be used to reverse the effect of a row/column-specific reserve()
// call. It removes all excessive capacity from the specified row (in case of a rowMajor matrix)
// or column (in case of a columnMajor matrix). The excessive capacity is assigned to the
// subsequent row/column.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::trim( size_t i ) const
{
   get().trim( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transposing the represented matrix.
//
// \return Reference to the transposed matrix.
//
// In case the access proxy represents a matrix-like data structure that provides an transpose()
// function, this function transposes the matrix.
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::transpose() const
{
   get().transpose();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the represented vector/matrix by the scalar value \a scalar.
//
// \param scalar The scalar value for the scaling.
// \return void
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// scale() function, this function performs a scaling by the given scalar.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the scalar value
inline void MatrixAccessProxy<MT>::scale( const Other& scalar ) const
{
   get().scale( scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  VECTOR/MATRIX LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific vector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// In case the access proxy represents a vector-like data structure that provides a find()
// function, this function can be used to check whether a specific element is contained in the
// sparse vector. It specifically searches for the element with index \a index. In case the
// element is found, the function returns an iterator to the element. Otherwise an iterator just
// past the last non-zero element of the compressed vector (the end() iterator) is returned. Note
// that the returned compressed vector iterator is subject to invalidation due to inserting
// operations via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator MatrixAccessProxy<MT>::find( size_t index ) const
{
   return get().find( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific matrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// In case the access proxy represents a matrix-like data structure that provides a find()
// function, this function can be used to check whether a specific element is contained in the
// sparse matrix. It specifically searches for the element with row index \a i and column index
// \a j. In case the element is found, the function returns an row/column iterator to the element.
// Otherwise an iterator just past the last non-zero element of row \a i or column \a j (the end()
// iterator) is returned. Note that the returned sparse matrix iterator is subject to invalidation
// due to inserting operations via the function call operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::find( size_t i, size_t j ) const
{
   return get().find( i, j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case the access proxy represents a vector-like data structure that provides a lowerBound()
// function, this function returns an iterator to the first element with an index not less then
// the given index. In combination with the upperBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned compressed
// vector iterator is subject to invalidation due to inserting operations via the subscript
// operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::lowerBound( size_t index ) const
{
   return get().lowerBound( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case the access proxy represents a matrix-like data structure that provides a lowerBound()
// function, this function can be used to search for a lower bound iterator within a row/column.
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index not less then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index not less then the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed matrix
// iterator is subject to invalidation due to inserting operations via the function call operator
// or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::lowerBound( size_t i, size_t j ) const
{
   return get().lowerBound( i, j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case the access proxy represents a vector-like data structure that provides a upperBound()
// function, this function returns an iterator to the first element with an index greater then
// the given index. In combination with the lowerBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned compressed
// vector iterator is subject to invalidation due to inserting operations via the subscript
// operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::upperBound( size_t index ) const
{
   return get().upperBound( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case the access proxy represents a matrix-like data structure that provides a upperBound()
// function, this function can be used to search for an upper bound iterator within a row/column.
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index greater then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index greater then the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed matrix
// iterator is subject to invalidation due to inserting operations via the function call operator
// or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Iterator
   MatrixAccessProxy<MT>::upperBound( size_t i, size_t j ) const
{
   return get().upperBound( i, j );
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
   typedef typename MatrixAccessProxy<MT1>::RawReference  LhsRawReference;
   typedef typename MatrixAccessProxy<MT2>::RawReference  RhsRawReference;
   return ( static_cast<LhsRawReference>( lhs ) == static_cast<RhsRawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( static_cast<RawReference>( lhs ) == rhs );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( lhs == static_cast<RawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT1>::RawReference  LhsRawReference;
   typedef typename MatrixAccessProxy<MT2>::RawReference  RhsRawReference;
   return ( static_cast<LhsRawReference>( lhs ) != static_cast<RhsRawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( static_cast<RawReference>( lhs ) != rhs );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( lhs != static_cast<RawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT1>::RawReference  LhsRawReference;
   typedef typename MatrixAccessProxy<MT2>::RawReference  RhsRawReference;
   return ( static_cast<LhsRawReference>( lhs ) < static_cast<RhsRawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( static_cast<RawReference>( lhs ) < rhs );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( lhs < static_cast<RawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT1>::RawReference  LhsRawReference;
   typedef typename MatrixAccessProxy<MT2>::RawReference  RhsRawReference;
   return ( static_cast<LhsRawReference>( lhs ) > static_cast<RhsRawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( static_cast<RawReference>( lhs ) > rhs );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( lhs > static_cast<RawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT1>::RawReference  LhsRawReference;
   typedef typename MatrixAccessProxy<MT2>::RawReference  RhsRawReference;
   return ( static_cast<LhsRawReference>( lhs ) <= static_cast<RhsRawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( static_cast<RawReference>( lhs ) <= rhs );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( lhs <= static_cast<RawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT1>::RawReference  LhsRawReference;
   typedef typename MatrixAccessProxy<MT2>::RawReference  RhsRawReference;
   return ( static_cast<LhsRawReference>( lhs ) >= static_cast<RhsRawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( static_cast<RawReference>( lhs ) >= rhs );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return ( lhs >= static_cast<RawReference>( rhs ) );
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
   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return os << static_cast<RawReference>( proxy );
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
inline typename MatrixAccessProxy<MT>::Iterator begin( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator begin( const MatrixAccessProxy<MT>& proxy, size_t i );

template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cbegin( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cbegin( const MatrixAccessProxy<MT>& proxy, size_t i );

template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator end( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator end( const MatrixAccessProxy<MT>& proxy, size_t i );

template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cend( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cend( const MatrixAccessProxy<MT>& proxy, size_t i );

template< typename MT >
inline size_t size( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline size_t rows( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline size_t columns( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline size_t capacity( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline size_t capacity( const MatrixAccessProxy<MT>& proxy, size_t i );

template< typename MT >
inline size_t nonZeros( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline size_t nonZeros( const MatrixAccessProxy<MT>& proxy, size_t i );

template< typename MT >
inline void reset( const MatrixAccessProxy<MT>& proxy );

template< typename MT >
inline void reset( const MatrixAccessProxy<MT>& proxy, size_t i );

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
/*!\brief Returns an iterator to the first element of the represented vector.
// \ingroup math
//
// \param proxy The given access proxy.
// \return Iterator to the first element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a begin()
// function, this function returns an iterator to the first element of the vector.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator begin( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a begin()
// function, this function returns an iterator to the first element of the of row/column \a i
// of the matrix. In case the given matrix is a row-major matrix the function returns an iterator
// to the first element of row \a i, in case it is a column-major matrix the function returns an
// iterator to the first element of column \a i.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator begin( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   return proxy.begin(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the represented vector.
// \ingroup math
//
// \param proxy The given access proxy.
// \return Iterator to the first element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a cbegin()
// function, this function returns an iterator to the first element of the vector.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cbegin( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.cbegin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a cbegin()
// function, this function returns an iterator to the first element of the of row/column \a i
// of the matrix. In case the given matrix is a row-major matrix the function returns an iterator
// to the first element of row \a i, in case it is a column-major matrix the function returns an
// iterator to the first element of column \a i.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cbegin( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   return proxy.cbegin(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the represented vector.
// \ingroup math
//
// \param proxy The given access proxy.
// \return Iterator just past the last element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides an end()
// function, this function returns an iterator just past the last element of the vector.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator end( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides an end()
// function, this function returns an iterator just past the last element of row/column \a i of
// the matrix. In case the given matrix is a row-major matrix the function returns an iterator
// just past the last element of row \a i, in case it is a column-major matrix the function
// returns an iterator just past the last element of column \a i.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::Iterator end( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   return proxy.end(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the represented vector.
// \ingroup math
//
// \param proxy The given access proxy.
// \return Iterator just past the last element of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a cend()
// function, this function returns an iterator just past the last element of the vector.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cend( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.cend();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a cend()
// function, this function returns an iterator just past the last element of row/column \a i of
// the matrix. In case the given matrix is a row-major matrix the function returns an iterator
// just past the last element of row \a i, in case it is a column-major matrix the function
// returns an iterator just past the last element of column \a i.
*/
template< typename MT >
inline typename MatrixAccessProxy<MT>::ConstIterator cend( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   return proxy.cend(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size/dimension of the represented vector.
// \ingroup math
//
// \param proxy The given access proxy.
// \return The size of the vector.
//
// In case the access proxy represents a vector-like data structure that provides a size()
// function, this function returns the current size/dimension of the vector.
*/
template< typename MT >
inline size_t size( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of rows of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \return The number of rows of the matrix.
//
// In case the access proxy represents a matrix-like data structure that provides a rows()
// function, this function returns the current number of rows of the matrix.
*/
template< typename MT >
inline size_t rows( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.rows();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \return The number of columns of the matrix.
//
// In case the access proxy represents a matrix-like data structure that provides a columns()
// function, this function returns the current number of columns of the matrix.
*/
template< typename MT >
inline size_t columns( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the represented vector/matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \return The capacity of the vector.
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// capacity() function, this function returns the current maximum capacity of the vector/matrix.
*/
template< typename MT >
inline size_t capacity( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.capacity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row/column of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a capacity()
// function, this function returns the current capacity of the specified row/column. In case
// the storage order is set to \a rowMajor the function returns the capacity of row \a i, in
// case the storage flag is set to \a columnMajor the function returns the capacity of column
// \a i.
*/
template< typename MT >
inline size_t capacity( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   return proxy.capacity(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the represented vector/matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \return The number of non-zero elements in the vector.
//
// In case the access proxy represents a vector- or matrix-like data structure that provides a
// nonZeros() function, this function returns the current number of non-zero elements in the
// vector/matrix. Note that the number of non-zero elements is always less than or equal to the
// current size of the vector.
*/
template< typename MT >
inline size_t nonZeros( const MatrixAccessProxy<MT>& proxy )
{
   return proxy.nonZeros();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// In case the access proxy represents a matrix-like data structure that provides a nonZeros()
// function, this function returns the current number of non-zero elements in the specified
// row/column. In case the storage order is set to \a rowMajor the function returns the number
// of non-zero elements in row \a i, in case the storage flag is set to \a columnMajor the
// function returns the number of non-zero elements in column \a i.
*/
template< typename MT >
inline size_t nonZeros( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   return proxy.nonZeros(i);
}
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
   proxy.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the represented matrix.
// \ingroup math
//
// \param proxy The given access proxy.
// \param i The index of the row/column to be resetted.
// \return void
//
// In case the access proxy represents a matrix-like data structure that provides a reset()
// function, this function resets all elements in the specified row/column of the given matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT >
inline void reset( const MatrixAccessProxy<MT>& proxy, size_t i )
{
   proxy.reset(i);
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
   proxy.clear();
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

   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   return isDefault( static_cast<RawReference>( proxy ) );
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

   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   swap( static_cast<RawReference>( a ), static_cast<RawReference>( b ) );
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

   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   swap( static_cast<RawReference>( a ), b );
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

   typedef typename MatrixAccessProxy<MT>::RawReference  RawReference;
   swap( a, static_cast<RawReference>( b ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
