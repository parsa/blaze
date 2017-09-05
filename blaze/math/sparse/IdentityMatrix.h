//=================================================================================================
/*!
//  \file blaze/math/sparse/IdentityMatrix.h
//  \brief Implementation of an identity matrix
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_SPARSE_IDENTITYMATRIX_H_
#define _BLAZE_MATH_SPARSE_IDENTITYMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DeclDiagTrait.h>
#include <blaze/math/traits/DeclHermTrait.h>
#include <blaze/math/traits/DeclLowTrait.h>
#include <blaze/math/traits/DeclSymTrait.h>
#include <blaze/math/traits/DeclUppTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/UnaryMapTrait.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup identity_matrix IdentityMatrix
// \ingroup sparse_matrix
*/
/*!\brief Efficient implementation of an \f$ N \times N \f$ identity matrix.
// \ingroup identity_matrix
//
// The IdentityMatrix class template is the representation of an immutable, arbitrary sized
// identity matrix with \f$ N \cdot N \f$ elements of arbitrary type. The type of the elements
// and the storage order of the matrix can be specified via the two template parameters:

   \code
   template< typename Type, bool SO >
   class IdentityMatrix;
   \endcode

//  - Type: specifies the type of the matrix elements. IdentityMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::rowMajor.
//
// It is not possible to insert, erase or modify the elements of an identity matrix. It is only
// possible to read from the elements:

   \code
   using blaze::rowMajor;

   // Creating a row-major 4x4 identity matrix with 4 rows and 4 columns
   IdentityMatrix<double,rowMajor> A( 4 );

   // The function call operator provides access to all possible elements of the identity matrix,
   // including the zero elements.
   A(1,2) = 2.0;       // Compilation error: It is not possible to write to an indentity matrix
   double d = A(2,1);  // Access to the element (2,1)

   // In order to traverse all non-zero elements currently stored in the matrix, the begin()
   // and end() functions can be used. In the example, all non-zero elements of the 2nd row
   // of A are traversed.
   for( IdentityMatrix<double,rowMajor>::Iterator i=A.begin(1); i!=A.end(1); ++i ) {
      ... = i->value();  // Access to the value of the non-zero element
      ... = i->index();  // Access to the index of the non-zero element
   }
   \endcode

// The use of IdentityMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, ...) can be performed on all possible combinations of row-major and column-major
// dense and sparse matrices with fitting element types. The following example gives an impression
// of the use of IdentityMatrix:

   \code
   using blaze::IdentityMatrix;
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   IdentityMatrix<double,rowMajor> A( 3 );  // Row-major 3x3 identity matrix

   DynamicMatrix<double,columnMajor> B( 3, 3 );  // Column-major 3x3 dynamic dense matrix
   CompressedMatrix<double,rowMajor> C( 3, 3 );  // Row-major 3x3 compressed sparse matrix
   CompressedMatrix<double,rowMajor> D( 3, 5 );  // Row-major 3x5 compressed sparse matrix
   // ... Initialization of B, C, and D

   DynamicMatrix<double,rowMajor>       E( A );  // Creation of a new row-major matrix as a copy of A
   CompressedMatrix<double,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;    // Addition of an identity matrix and a dense matrix
   E = C - A;    // Subtraction of a sparse matrix and an identity matrix
   F = A * D;    // Matrix multiplication between two matrices of different element types

   E = 2.0 * A;  // Scaling of an identity matrix
   F = A * 2.0;  // Scaling of an identity matrix
   \endcode
*/
template< typename Type                    // Data type of the matrix
        , bool SO = defaultStorageOrder >  // Storage order
class IdentityMatrix
   : public SparseMatrix< IdentityMatrix<Type,SO>, SO >
{
 public:
   //**Type definitions****************************************************************************
   using This           = IdentityMatrix<Type,SO>;   //!< Type of this IdentityMatrix instance.
   using BaseType       = SparseMatrix<This,SO>;     //!< Base type of this IdentityMatrix instance.
   using ResultType     = This;                      //!< Result type for expression template evaluations.
   using OppositeType   = IdentityMatrix<Type,!SO>;  //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType  = IdentityMatrix<Type,!SO>;  //!< Transpose type for expression template evaluations.
   using ElementType    = Type;                      //!< Type of the identity matrix elements.
   using ReturnType     = const Type;                //!< Return type for expression template evaluations.
   using CompositeType  = const This&;               //!< Data type for composite expression templates.
   using Reference      = const Type;                //!< Reference to a identity matrix element.
   using ConstReference = const Type;                //!< Reference to a constant identity matrix element.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain an IdentityMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using Other = IdentityMatrix<NewType,SO>;  //!< The type of the other IdentityMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a IdentityMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using Other = IdentityMatrix<Type,SO>;  //!< The type of the other IdentityMatrix.
   };
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the identity matrix.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the identity matrix.
      using Element = ValueIndexPair<Type>;

      using IteratorCategory = std::forward_iterator_tag;  //!< The iterator category.
      using ValueType        = Element;                    //!< Type of the underlying pointers.
      using PointerType      = ValueType*;                 //!< Pointer return type.
      using ReferenceType    = ValueType&;                 //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                  //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying pointers.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Default constructor**********************************************************************
      /*!\brief Default constructor for the ConstIterator class.
      */
      inline ConstIterator()
         : index_()  // Index to the current identity matrix element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param index Index to the initial matrix element.
      */
      inline ConstIterator( size_t index )
         : index_( index )  // Index to the current identity matrix element
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ConstIterator& operator++() {
         ++index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline ConstIterator operator++( int ) {
         ConstIterator tmp( *this );
         ++index_;
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline const Element operator*() const {
         return Element( Type(1), index_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return Reference to the sparse matrix element at the current iterator position.
      */
      inline const ConstIterator* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse element.
      //
      // \return The current value of the sparse element.
      */
      inline Type value() const {
         return Type(1);
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return index_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side ConstIterator object.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side ConstIterator object.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two ConstIterator objects.
      //
      // \param rhs The right-hand side ConstIterator object.
      // \return The number of elements between the two ConstIterator objects.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return index_ - rhs.index_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      size_t index_;  //!< Index to the current identity matrix element.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using Iterator = ConstIterator;  //!< Iterator over non-constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   enum : bool { smpAssignable = !IsSMPAssignable<Type>::value };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline IdentityMatrix() noexcept;
   explicit inline IdentityMatrix( size_t n ) noexcept;

   template< typename MT, bool SO2 >
   explicit inline IdentityMatrix( const Matrix<MT,SO2>& m );

   // No explicitly declared copy constructor.
   // No explicitly declared move constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline ConstReference at( size_t i, size_t j ) const;
   inline ConstIterator  begin ( size_t i ) const noexcept;
   inline ConstIterator  cbegin( size_t i ) const noexcept;
   inline ConstIterator  end   ( size_t i ) const noexcept;
   inline ConstIterator  cend  ( size_t i ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   template< typename MT, bool SO2 >
   inline IdentityMatrix& operator=( const Matrix<MT,SO2>& rhs );

   // No explicitly declared copy assignment operator.
   // No explicitly declared move assignment operator.
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t rows() const noexcept;
   inline size_t columns() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   clear();
          void   resize( size_t n );
   inline void   swap( IdentityMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline IdentityMatrix& transpose();
   inline IdentityMatrix& ctranspose();
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t n_;  //!< The current number of rows and columns of the identity matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE      ( Type );
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
/*!\brief The default constructor for IdentityMatrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline IdentityMatrix<Type,SO>::IdentityMatrix() noexcept
   : n_( 0UL )  // The current number of rows and columns of the identity matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for an identity matrix of size \f$ N \times N \f$.
//
// \param n The number of rows and columns of the matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline IdentityMatrix<Type,SO>::IdentityMatrix( size_t n ) noexcept
   : n_( n )  // The current number of rows and columns of the identity matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor for different identity matrices.
//
// \param m Identity matrix to be copied.
// \exception std::invalid_argument Invalid setup of identity matrix.
//
// The matrix is sized according to the given \f$ N \times N \f$ identity matrix and
// initialized as a copy of this matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
template< typename MT    // Type of the foreign identity matrix
        , bool SO2 >     // Storage order of the foreign identity matrix
inline IdentityMatrix<Type,SO>::IdentityMatrix( const Matrix<MT,SO2>& m )
   : n_( (~m).rows() )  // The current number of rows and columns of the identity matrix
{
   if( !IsIdentity<MT>::value && !isIdentity( ~m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of identity matrix" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the identity matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstReference
   IdentityMatrix<Type,SO>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid identity matrix column access index" );

   if( i == j )
      return Type( 1 );
   else
      return Type( 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstReference
   IdentityMatrix<Type,SO>::at( size_t i, size_t j ) const
{
   if( i >= n_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= n_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::begin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::cbegin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::end( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i+1UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::cend( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i+1UL );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Assignment operator for different identity matrices.
//
// \param rhs Identity matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid assignment to identity matrix.
//
// The matrix is resized according to the given \f$ N \times N \f$ identity matrix and
// initialized as a copy of this matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side identity matrix
        , bool SO2 >     // Storage order of the right-hand side identity matrix
inline IdentityMatrix<Type,SO>&
   IdentityMatrix<Type,SO>::operator=( const Matrix<MT,SO2>& rhs )
{
   if( !IsIdentity<MT>::value && !isIdentity( ~rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment of identity matrix" );
   }

   n_ = (~rhs).rows();

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the identity matrix.
//
// \return The number of rows of the identity matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t IdentityMatrix<Type,SO>::rows() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the identity matrix.
//
// \return The number of columns of the identity matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t IdentityMatrix<Type,SO>::columns() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the identity matrix.
//
// \return The capacity of the identity matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t IdentityMatrix<Type,SO>::capacity() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row/column.
//
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// This function returns the current capacity of the specified row/column. In case the
// storage order is set to \a rowMajor the function returns the capacity of row \a i,
// in case the storage flag is set to \a columnMajor the function returns the capacity
// of column \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t IdentityMatrix<Type,SO>::capacity( size_t i ) const noexcept
{
   UNUSED_PARAMETER( i );

   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return 1UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the identity matrix
//
// \return The number of non-zero elements in the identity matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t IdentityMatrix<Type,SO>::nonZeros() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column.
// In case the storage order is set to \a rowMajor the function returns the number of non-zero
// elements in row \a i, in case the storage flag is set to \a columnMajor the function returns
// the number of non-zero elements in column \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t IdentityMatrix<Type,SO>::nonZeros( size_t i ) const
{
   UNUSED_PARAMETER( i );

   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return 1UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the identity matrix.
//
// \return void
//
// After the clear() function, the size of the identity matrix is 0.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void IdentityMatrix<Type,SO>::clear()
{
   n_ = 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the identity matrix.
//
// \param n The new number of rows and columns of the identity matrix.
// \return void
//
// This function resizes the matrix using the given size to \f$ n \times n \f$. Note that this
// function may invalidate all existing views (submatrices, rows, columns, ...) on the matrix if
// it is used to shrink the matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
void IdentityMatrix<Type,SO>::resize( size_t n )
{
   n_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two sparse matrices.
//
// \param m The identity matrix to be swapped.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void IdentityMatrix<Type,SO>::swap( IdentityMatrix& m ) noexcept
{
   std::swap( n_, m.n_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific matrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// matrix. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an row/column iterator to the element.
// Otherwise an iterator just past the last non-zero element of row \a i or column \a j (the
// end() iterator) is returned.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::find( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( SO  || i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( !SO || j < columns(), "Invalid identity matrix column access index" );

   if( i == j )
      return begin( i );
   else
      return end( SO ? j : i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index not less then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index not less then the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::lowerBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( SO  || i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( !SO || j < columns(), "Invalid identity matrix column access index" );

   if( ( !SO && j <= i ) || ( SO && i <= j ) )
      return begin( SO ? j : i );
   else
      return end( SO ? j : i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index greater then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index greater then the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename IdentityMatrix<Type,SO>::ConstIterator
   IdentityMatrix<Type,SO>::upperBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( SO  || i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( !SO || j < columns(), "Invalid identity matrix column access index" );

   if( ( !SO && j < i ) || ( SO && i < j ) )
      return begin( SO ? j : i );
   else
      return end( SO ? j : i );
}
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief In-place transpose of the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline IdentityMatrix<Type,SO>& IdentityMatrix<Type,SO>::transpose()
{
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place conjugate transpose of the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline IdentityMatrix<Type,SO>& IdentityMatrix<Type,SO>::ctranspose()
{
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the matrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address can alias with the matrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool IdentityMatrix<Type,SO>::canAlias( const Other* alias ) const noexcept
{
   UNUSED_PARAMETER( alias );

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address is aliased with the matrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool IdentityMatrix<Type,SO>::isAliased( const Other* alias ) const noexcept
{
   UNUSED_PARAMETER( alias );

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix can be used in SMP assignments.
//
// \return \a true in case the matrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the matrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline bool IdentityMatrix<Type,SO>::canSMPAssign() const noexcept
{
   return false;
}
//*************************************************************************************************








//=================================================================================================
//
//  IDENTITYMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name IdentityMatrix operators */
//@{
template< typename Type, bool SO >
inline void reset( IdentityMatrix<Type,SO>& m );

template< typename Type, bool SO >
inline void reset( IdentityMatrix<Type,SO>& m, size_t i );

template< typename Type, bool SO >
inline void clear( IdentityMatrix<Type,SO>& m );

template< bool RF, typename Type, bool SO >
inline bool isDefault( const IdentityMatrix<Type,SO>& m );

template< typename Type, bool SO >
inline bool isIntact( const IdentityMatrix<Type,SO>& m );

template< typename Type, bool SO >
inline void swap( IdentityMatrix<Type,SO>& a, IdentityMatrix<Type,SO>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given identity matrix.
// \ingroup identity_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void reset( IdentityMatrix<Type,SO>& m )
{
   UNUSED_PARAMETER( m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given identity matrix.
// \ingroup identity_matrix
//
// \param m The matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given identity matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void reset( IdentityMatrix<Type,SO>& m, size_t i )
{
   UNUSED_PARAMETER( m, i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given identity matrix.
// \ingroup identity_matrix
//
// \param m The matrix to be cleared.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void clear( IdentityMatrix<Type,SO>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given identity matrix is in default state.
// \ingroup identity_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix's rows and columns are zero, \a false otherwise.
//
// This function checks whether the identity matrix is in default (constructed) state, i.e. if
// it's number of rows and columns is 0. In case it is in default state, the function returns
// \a true, else it will return \a false. The following example demonstrates the use of the
// \a isDefault() function:

   \code
   blaze::IdentityMatrix<int> I;
   // ... Resizing and initialization
   if( isDefault( I ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( I ) ) { ... }
   \endcode
*/
template< bool RF        // Relaxation flag
        , typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline bool isDefault( const IdentityMatrix<Type,SO>& m )
{
   return ( m.rows() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given identity matrix are intact.
// \ingroup identity_matrix
//
// \param m The identity matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the identity matrix are intact, i.e. if
// its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::IdentityMatrix<int> I;
   // ... Resizing and initialization
   if( isIntact( I ) ) { ... }
   \endcode
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline bool isIntact( const IdentityMatrix<Type,SO>& m )
{
   UNUSED_PARAMETER( m );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two identity matrices.
// \ingroup identity_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void swap( IdentityMatrix<Type,SO>& a, IdentityMatrix<Type,SO>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of an identity matrix and a dense vector
//        (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup identity_matrix
//
// \param mat The left-hand side identity matrix for the multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a row-major sparse matrix and a dense
// vector:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::IdentityMatrix<double,rowMajor> A;
   blaze::DynamicVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns a reference to the given dense vector. In case the current size of
// the vector \a vec doesn't match the current number of columns of the matrix \a mat, a
// \a std::invalid_argument is thrown.
*/
template< typename T   // Data type of the left-hand side identity matrix
        , bool SO      // Storage order of the left-hand side identity matrix
        , typename VT  // Type of the right-hand side dense vector
        , typename = EnableIf_< IsSame< T, ElementType_<VT> > > >
inline decltype(auto)
   operator*( const IdentityMatrix<T,SO>& mat, const DenseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   if( (~mat).columns() != (~vec).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix and vector sizes do not match" );
   }

   return (~vec);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose dense vector and an
//        identity matrix (\f$ \vec{y}^T=\vec{x}^T*A \f$).
// \ingroup identity_matrix
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side identity matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose dense vector and an identity
// matrix:

   \code
   using blaze::rowVector;
   using blaze::rowMajor;

   blaze::DynamicVector<double,rowVector> x, y;
   blaze::IentityMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// The operator returns a reference to the given dense vector. In case the current size of
// the vector \a vec doesn't match the current number of rows of the matrix \a mat, a
// \a std::invalid_argument is thrown.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename T   // Data type of the right-hand side identity matrix
        , bool SO      // Storage order of the right-hand side identity matrix
        , typename = EnableIf_< IsSame< ElementType_<VT>, T > > >
inline decltype(auto)
   operator*( const DenseVector<VT,true>& vec, const IdentityMatrix<T,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   if( (~vec).size() != (~mat).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector and matrix sizes do not match" );
   }

   return (~vec);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a identity matrix and a sparse vector
//        (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup identity_matrix
//
// \param mat The left-hand side identity matrix for the multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a row-major sparse matrix and a sparse
// vector:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::IdentityMatrix<double,rowMajor> A;
   blaze::CompressedVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns a reference to the given sparse vector. In case the current size of
// the vector \a vec doesn't match the current number of columns of the matrix \a mat, a
// \a std::invalid_argument is thrown.
*/
template< typename T   // Data type of the left-hand side identity matrix
        , bool SO      // Storage order of the left-hand side identity matrix
        , typename VT  // Type of the right-hand side sparse vector
        , typename = EnableIf_< IsSame< T, ElementType_<VT> > > >
inline decltype(auto)
   operator*( const IdentityMatrix<T,SO>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   if( (~mat).columns() != (~vec).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix and vector sizes do not match" );
   }

   return (~vec);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose sparse vector and an
//        identity matrix (\f$ \vec{y}^T=\vec{x}^T*A \f$).
// \ingroup identity_matrix
//
// \param vec The left-hand side transpose sparse vector for the multiplication.
// \param mat The right-hand side identity matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose sparse vector and an identity
// matrix:

   \code
   using blaze::rowVector;
   using blaze::rowMajor;

   blaze::CompressedVector<double,rowVector> x, y;
   blaze::IentityMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// The operator returns a reference to the given sparse vector. In case the current size of
// the vector \a vec doesn't match the current number of rows of the matrix \a mat, a
// \a std::invalid_argument is thrown.
*/
template< typename VT  // Type of the left-hand side sparse vector
        , typename T   // Data type of the right-hand side identity matrix
        , bool SO      // Storage order of the right-hand side identity matrix
        , typename = EnableIf_< IsSame< ElementType_<VT>, T > > >
inline decltype(auto)
   operator*( const SparseVector<VT,true>& vec, const IdentityMatrix<T,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   if( (~vec).size() != (~mat).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector and matrix sizes do not match" );
   }

   return (~vec);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of an identity matrix and a dense
//        matrix (\f$ A=B*C \f$).
// \ingroup identity_matrix
//
// \param lhs The left-hand side identity matrix for the multiplication.
// \param rhs The right-hand side dense matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of an identity matrix and a dense matrix:

   \code
   using blaze::rowMajor;

   blaze::IdentityMatrix<double,rowMajor> A;
   blaze::DynamicMatrix<double,rowMajor> B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns a reference to the given dense matrix. In case the current sizes of the
// two given matrices don't match, a \a std::invalid_argument is thrown.
*/
template< typename T   // Data type of the left-hand side identity matrix
        , bool SO1     // Storage order of the left-hand side identity matrix
        , typename MT  // Type of the right-hand side dense matrix
        , bool SO2     // Storage order of the right-hand side dense matrix
        , typename = EnableIf_< IsSame< T, ElementType_<MT> > > >
inline decltype(auto)
   operator*( const IdentityMatrix<T,SO1>& lhs, const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return (~rhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix and an identity matrix
//        (\f$ A=B*C \f$).
// \ingroup identity_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side identity matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a dense matrix and an identity matrix:

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> A, C;
   blaze::IdentityMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns a reference to the given dense matrix. In case the current sizes of the
// two given matrices don't match, a \a std::invalid_argument is thrown.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , bool SO1     // Storage order of the left-hand side dense matrix
        , typename T   // Data type of the right-hand side identity matrix
        , bool SO2     // Storage order of the right-hand side identity matrix
        , typename = EnableIf_< IsSame< ElementType_<MT>, T > > >
inline decltype(auto)
   operator*( const DenseMatrix<MT,SO1>& lhs, const IdentityMatrix<T,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return (~lhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of an identity matrix and a sparse
//        matrix (\f$ A=B*C \f$).
// \ingroup identity_matrix
//
// \param lhs The left-hand side identity matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of an identity matrix and a sparse matrix:

   \code
   using blaze::rowMajor;

   blaze::IdentityMatrix<double,rowMajor> A;
   blaze::DynamicMatrix<double,rowMajor> B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns a reference to the given sparse matrix. In case the current sizes of the
// two given matrices don't match, a \a std::invalid_argument is thrown.
*/
template< typename T   // Data type of the left-hand side identity matrix
        , bool SO1     // Storage order of the left-hand side identity matrix
        , typename MT  // Type of the right-hand side sparse matrix
        , bool SO2     // Storage order of the right-hand side sparse matrix
        , typename = EnableIf_< IsSame< T, ElementType_<MT> > > >
inline decltype(auto)
   operator*( const IdentityMatrix<T,SO1>& lhs, const SparseMatrix<MT,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return (~rhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse matrix and an identity matrix
//        (\f$ A=B*C \f$).
// \ingroup identity_matrix
//
// \param lhs The left-hand side sparse matrix for the multiplication.
// \param rhs The right-hand side identity matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a sparse matrix and an identity matrix:

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> A, C;
   blaze::IdentityMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns a reference to the given sparse matrix. In case the current sizes of the
// two given matrices don't match, a \a std::invalid_argument is thrown.
*/
template< typename MT  // Type of the left-hand side sparse matrix
        , bool SO1     // Storage order of the left-hand side sparse matrix
        , typename T   // Data type of the right-hand side identity matrix
        , bool SO2     // Storage order of the right-hand side identity matrix
        , typename = EnableIf_< IsSame< ElementType_<MT>, T > > >
inline decltype(auto)
   operator*( const SparseMatrix<MT,SO1>& lhs, const IdentityMatrix<T,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return (~lhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of two identity matrices (\f$ A=B*C \f$).
// \ingroup identity_matrix
//
// \param lhs The left-hand side identity matrix for the multiplication.
// \param rhs The right-hand side identity matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of two identity matrices:

   \code
   using blaze::rowMajor;

   blaze::IdentityMatrix<double,rowMajor> A, B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns an identity matrix. In case the current sizes of the two given matrices
// don't match, a \a std::invalid_argument is thrown.
*/
template< typename T1  // Data type of the left-hand side identity matrix
        , bool SO1     // Storage order of the left-hand side identity matrix
        , typename T2  // Data type of the right-hand side dense matrix
        , bool SO2 >   // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator*( const IdentityMatrix<T1,SO1>& lhs, const IdentityMatrix<T2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return IdentityMatrix< MultTrait_<T1,T2>, SO1 >( (~lhs).rows() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Declares the given matrix expression \a m as identity matrix.
// \ingroup identity_matrix
//
// \param m The input matrix.
// \return The redeclared matrix.
// \exception std::invalid_argument Invalid identity matrix specification.
//
// The \a declid function declares the given dense or sparse matrix expression \a m as identity
// matrix. In case the given matrix is not a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a declid function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = declid( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline IdentityMatrix<ElementType_<MT>,SO>
   declid( const Matrix<MT,SO>& m )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( ~m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid identity matrix specification" );
   }

   return IdentityMatrix<ElementType_<MT>,SO>( (~m).rows() );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISSQUARE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsSquare< IdentityMatrix<MT,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsSymmetric< IdentityMatrix<MT,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsHermitian< IdentityMatrix<MT,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsUniLower< IdentityMatrix<MT,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsUniUpper< IdentityMatrix<MT,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct IsResizable< IdentityMatrix<T,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2, size_t M, size_t N >
struct AddTrait< IdentityMatrix<T1,SO>, StaticMatrix<T2,M,N,SO> >
{
   using Type = StaticMatrix< AddTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct AddTrait< IdentityMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   using Type = StaticMatrix< AddTrait_<T1,T2>, M, N, SO2 >;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct AddTrait< StaticMatrix<T1,M,N,SO>, IdentityMatrix<T2,SO> >
{
   using Type = StaticMatrix< AddTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct AddTrait< StaticMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = StaticMatrix< AddTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, bool SO, typename T2, size_t M, size_t N >
struct AddTrait< IdentityMatrix<T1,SO>, HybridMatrix<T2,M,N,SO> >
{
   using Type = HybridMatrix< AddTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct AddTrait< IdentityMatrix<T1,SO1>, HybridMatrix<T2,M,N,SO2> >
{
   using Type = HybridMatrix< AddTrait_<T1,T2>, M, N, SO2 >;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct AddTrait< HybridMatrix<T1,M,N,SO>, IdentityMatrix<T2,SO> >
{
   using Type = HybridMatrix< AddTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct AddTrait< HybridMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = HybridMatrix< AddTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< IdentityMatrix<T1,SO>, DynamicMatrix<T2,SO> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< IdentityMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< DynamicMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< DynamicMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO, typename T2, bool AF, bool PF >
struct AddTrait< IdentityMatrix<T1,SO>, CustomMatrix<T2,AF,PF,SO> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool AF, bool PF, bool SO2 >
struct AddTrait< IdentityMatrix<T1,SO1>, CustomMatrix<T2,AF,PF,SO2> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool AF, bool PF, bool SO, typename T2 >
struct AddTrait< CustomMatrix<T1,AF,PF,SO>, IdentityMatrix<T2,SO> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool AF, bool PF, bool SO1, typename T2, bool SO2 >
struct AddTrait< CustomMatrix<T1,AF,PF,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< AddTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< IdentityMatrix<T1,SO>, CompressedMatrix<T2,SO> >
{
   using Type = CompressedMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< IdentityMatrix<T1,SO1>, CompressedMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< AddTrait_<T1,T2>, false >;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< CompressedMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = CompressedMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< CompressedMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< AddTrait_<T1,T2>, false >;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< IdentityMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = CompressedMatrix< AddTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< IdentityMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< AddTrait_<T1,T2>, false >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2, size_t M, size_t N >
struct SubTrait< IdentityMatrix<T1,SO>, StaticMatrix<T2,M,N,SO> >
{
   using Type = StaticMatrix< SubTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct SubTrait< IdentityMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   using Type = StaticMatrix< SubTrait_<T1,T2>, M, N, SO2 >;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct SubTrait< StaticMatrix<T1,M,N,SO>, IdentityMatrix<T2,SO> >
{
   using Type = StaticMatrix< SubTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct SubTrait< StaticMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = StaticMatrix< SubTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, bool SO, typename T2, size_t M, size_t N >
struct SubTrait< IdentityMatrix<T1,SO>, HybridMatrix<T2,M,N,SO> >
{
   using Type = HybridMatrix< SubTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct SubTrait< IdentityMatrix<T1,SO1>, HybridMatrix<T2,M,N,SO2> >
{
   using Type = HybridMatrix< SubTrait_<T1,T2>, M, N, SO2 >;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct SubTrait< HybridMatrix<T1,M,N,SO>, IdentityMatrix<T2,SO> >
{
   using Type = HybridMatrix< SubTrait_<T1,T2>, M, N, SO >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct SubTrait< HybridMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = HybridMatrix< SubTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< IdentityMatrix<T1,SO>, DynamicMatrix<T2,SO> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< IdentityMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< DynamicMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< DynamicMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO, typename T2, bool AF, bool PF >
struct SubTrait< IdentityMatrix<T1,SO>, CustomMatrix<T2,AF,PF,SO> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool AF, bool PF, bool SO2 >
struct SubTrait< IdentityMatrix<T1,SO1>, CustomMatrix<T2,AF,PF,SO2> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool AF, bool PF, bool SO, typename T2 >
struct SubTrait< CustomMatrix<T1,AF,PF,SO>, IdentityMatrix<T2,SO> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO >;
};

template< typename T1, bool AF, bool PF, bool SO1, typename T2, bool SO2 >
struct SubTrait< CustomMatrix<T1,AF,PF,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< SubTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< IdentityMatrix<T1,SO>, CompressedMatrix<T2,SO> >
{
   using Type = CompressedMatrix< SubTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< IdentityMatrix<T1,SO1>, CompressedMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< SubTrait_<T1,T2>, false >;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< CompressedMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = CompressedMatrix< SubTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< CompressedMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< SubTrait_<T1,T2>, false >;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< IdentityMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = CompressedMatrix< SubTrait_<T1,T2> , SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< IdentityMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< SubTrait_<T1,T2> , false >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHURTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct SchurTrait< IdentityMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct SchurTrait< StaticMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct SchurTrait< IdentityMatrix<T1,SO1>, HybridMatrix<T2,M,N,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct SchurTrait< HybridMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SchurTrait< IdentityMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SchurTrait< DynamicMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool SO1, typename T2, bool AF, bool PF, bool SO2 >
struct SchurTrait< IdentityMatrix<T1,SO1>, CustomMatrix<T2,AF,PF,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool AF, bool PF, bool SO1, typename T2, bool SO2 >
struct SchurTrait< CustomMatrix<T1,AF,PF,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO2 >;
};

template< typename T1, bool SO, typename T2 >
struct SchurTrait< IdentityMatrix<T1,SO>, CompressedMatrix<T2,SO> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SchurTrait< IdentityMatrix<T1,SO1>, CompressedMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, false >;
};

template< typename T1, bool SO, typename T2 >
struct SchurTrait< CompressedMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SchurTrait< CompressedMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, false >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SchurTrait< IdentityMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = IdentityMatrix< MultTrait_<T1,T2>, SO1 >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2 >
struct MultTrait< IdentityMatrix<T1,SO>, T2, EnableIf_< IsNumeric<T2> > >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO >;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< T1, IdentityMatrix<T2,SO>, EnableIf_< IsNumeric<T1> > >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO >;
};

template< typename T1, bool SO, typename T2, size_t N >
struct MultTrait< IdentityMatrix<T1,SO>, StaticVector<T2,N,false> >
{
   using Type = StaticVector< MultTrait_<T1,T2>, N, false >;
};

template< typename T1, size_t N, typename T2, bool SO >
struct MultTrait< StaticVector<T1,N,true>, IdentityMatrix<T2,SO> >
{
   using Type = StaticVector< MultTrait_<T1,T2>, N, true >;
};

template< typename T1, bool SO, typename T2, size_t N >
struct MultTrait< IdentityMatrix<T1,SO>, HybridVector<T2,N,false> >
{
   using Type = HybridVector< MultTrait_<T1,T2>, N, false >;
};

template< typename T1, size_t N, typename T2, bool SO >
struct MultTrait< HybridVector<T1,N,true>, IdentityMatrix<T2,SO> >
{
   using Type = HybridVector< MultTrait_<T1,T2>, N, true >;
};

template< typename T1, bool SO, typename T2 >
struct MultTrait< IdentityMatrix<T1,SO>, DynamicVector<T2,false> >
{
   using Type = DynamicVector< MultTrait_<T1,T2>, false >;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< DynamicVector<T1,true>, IdentityMatrix<T2,SO> >
{
   using Type = DynamicVector< MultTrait_<T1,T2>, true >;
};

template< typename T1, bool SO, typename T2, bool AF, bool PF >
struct MultTrait< IdentityMatrix<T1,SO>, CustomVector<T2,AF,PF,false> >
{
   using Type = DynamicVector< MultTrait_<T1,T2>, false >;
};

template< typename T1, bool AF, bool PF, typename T2, bool SO >
struct MultTrait< CustomVector<T1,AF,PF,true>, IdentityMatrix<T2,SO> >
{
   using Type = DynamicVector< MultTrait_<T1,T2>, true >;
};

template< typename T1, bool SO, typename T2 >
struct MultTrait< IdentityMatrix<T1,SO>, CompressedVector<T2,false> >
{
   using Type = CompressedVector< MultTrait_<T1,T2>, false >;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< CompressedVector<T1,true>, IdentityMatrix<T2,SO> >
{
   using Type = CompressedVector< MultTrait_<T1,T2>, true >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct MultTrait< IdentityMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   using Type = StaticMatrix< MultTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct MultTrait< StaticMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = StaticMatrix< MultTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct MultTrait< IdentityMatrix<T1,SO1>, HybridMatrix<T2,M,N,SO2> >
{
   using Type = HybridMatrix< MultTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct MultTrait< HybridMatrix<T1,M,N,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = HybridMatrix< MultTrait_<T1,T2>, M, N, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< IdentityMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< DynamicMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool AF, bool PF, bool SO2 >
struct MultTrait< IdentityMatrix<T1,SO1>, CustomMatrix<T2,AF,PF,SO2> >
{
   using Type = DynamicMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool AF, bool PF, bool SO1, typename T2, bool SO2 >
struct MultTrait< CustomMatrix<T1,AF,PF,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = DynamicMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< IdentityMatrix<T1,SO1>, CompressedMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< CompressedMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = CompressedMatrix< MultTrait_<T1,T2>, SO1 >;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< IdentityMatrix<T1,SO1>, IdentityMatrix<T2,SO2> >
{
   using Type = IdentityMatrix< MultTrait_<T1,T2>, SO1 >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2 >
struct DivTrait< IdentityMatrix<T1,SO>, T2, EnableIf_< IsNumeric<T2> > >
{
   using Type = CompressedMatrix< DivTrait_<T1,T2>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UNARYMAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename OP >
struct UnaryMapTrait< IdentityMatrix<T,SO>, OP >
{
   using Type = CompressedMatrix< UnaryMapTrait_<T,OP>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Abs >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Abs>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Floor >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Floor>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Ceil >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Ceil>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Trunc >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Trunc>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Round >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Round>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Conj >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Conj>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Real >
{
   using Type = IdentityMatrix< UnaryMapTrait_<T,Real>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Sqrt >
{
   using Type = IdentityMatrix< UnaryMapTrait<T,Sqrt>, SO >;
};

template< typename T, bool SO >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Cbrt >
{
   using Type = IdentityMatrix< UnaryMapTrait<T,Cbrt>, SO >;
};

template< typename T, bool SO, typename ET >
struct UnaryMapTrait< IdentityMatrix<T,SO>, Pow<ET> >
{
   using Type = IdentityMatrix< UnaryMapTrait< T, Pow<ET> >, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSYMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct DeclSymTrait< IdentityMatrix<T,SO> >
{
   using Type = IdentityMatrix<T,SO>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLHERMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct DeclHermTrait< IdentityMatrix<T,SO> >
{
   using Type = IdentityMatrix<T,SO>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct DeclLowTrait< IdentityMatrix<T,SO> >
{
   using Type = IdentityMatrix<T,SO>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct DeclUppTrait< IdentityMatrix<T,SO> >
{
   using Type = IdentityMatrix<T,SO>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLDIAGTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct DeclDiagTrait< IdentityMatrix<T,SO> >
{
   using Type = IdentityMatrix<T,SO>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2 >
struct HighType< IdentityMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = IdentityMatrix< typename HighType<T1,T2>::Type, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2 >
struct LowType< IdentityMatrix<T1,SO>, IdentityMatrix<T2,SO> >
{
   using Type = IdentityMatrix< typename LowType<T1,T2>::Type, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct SubmatrixTrait< IdentityMatrix<T,SO> >
{
   using Type = CompressedMatrix<T,SO>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, size_t... RAs >
struct RowTrait< IdentityMatrix<T,SO>, RAs... >
{
   using Type = CompressedVector<T,true>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, size_t... CAs >
struct ColumnTrait< IdentityMatrix<T,SO>, CAs... >
{
   using Type = CompressedVector<T,false>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, ptrdiff_t... BAs >
struct BandTrait< IdentityMatrix<T,SO>, BAs... >
{
   using Type = CompressedVector<T,defaultTransposeFlag>;
};

template< typename T, bool SO >
struct BandTrait< IdentityMatrix<T,SO>, 0L >
{
   using Type = DynamicVector<T,defaultTransposeFlag>;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
