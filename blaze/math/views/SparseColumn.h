//=================================================================================================
/*!
//  \file blaze/math/views/SparseColumn.h
//  \brief Header file for the SparseColumn class template
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

#ifndef _BLAZE_MATH_VIEWS_SPARSECOLUMN_H_
#define _BLAZE_MATH_VIEWS_SPARSECOLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/expressions/Column.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Functions.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Null.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <blaze/util/Unused.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup sparse_column SparseColumn
// \ingroup views
*/
/*!\brief Reference to a specific column of a sparse matrix.
// \ingroup sparse_column
//
// The SparseColumn template represents a reference to a specific column of a sparse matrix
// primitive. The type of the sparse matrix is specified via the first template parameter:

   \code
   template< typename MT, bool SO >
   class SparseColumn;
   \endcode

//  - MT: specifies the type of the sparse matrix primitive. SparseColumn can be used with every
//        sparse matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the sparse matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//  - SF: specifies whether the given matrix is a symmetric matrix or not. Also this parameter
//        doesn't have to be explicitly defined, but is automatically derived from the first
//        template parameter.
//
//
// \n \section sparse_column_setup Setup of Sparse Columns
//
// A reference to a sparse column can be created very conveniently via the \c column() function.
// This reference can be treated as any other column vector, i.e. it can be assigned to, it can
// be copied from, and it can be used in arithmetic operations. The reference can also be used on
// both sides of an assignment: The column can either be used as an alias to grant write access
// to a specific column of a matrix primitive on the left-hand side of an assignment or to grant
// read-access to a specific column of a matrix primitive or expression on the right-hand side
// of an assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicVector<double,blaze::columnVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::columnVector>  SparseVectorType;
   typedef blaze::CompressedMatrix<double,blaze::columnMajor>   SparseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   SparseMatrixType A, B;
   // ... Resizing and initialization

   // Setting the 2nd column of matrix A to x
   blaze::SparseColumn<SparseMatrixType> col2 = column( A, 2UL );
   col2 = x;

   // Setting the 3rd column of matrix B to y
   column( B, 3UL ) = y;

   // Setting x to the 1st column of matrix B
   x = column( B, 1UL );

   // Setting y to the 4th column of the result of the matrix multiplication
   y = column( A * B, 4UL );
   \endcode

// \n \section sparse_column_element_access Element access
//
// A sparse column can be used like any other column vector. For instance, the elements of the
// sparse column can be directly accessed with the subscript operator.

   \code
   typedef blaze::CompressedMatrix<double,blaze::columnMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the 4th column of matrix A
   blaze::SparseColumn<MatrixType> col4 = column( A, 4UL );

   // Setting the 1st element of the sparse column, which corresponds
   // to the 1st element in the 4th column of matrix A
   col4[1] = 2.0;
   \endcode

// The numbering of the column elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of rows of the referenced matrix. Alternatively, the elements of a
// column can be traversed via iterators. Just as with vectors, in case of non-const columns,
// \c begin() and \c end() return an Iterator, which allows a manipulation of the non-zero
// values, in case of constant columns a ConstIterator is returned:

   \code
   typedef blaze::CompressedMatrix<int,blaze::columnMajor>  MatrixType;
   typedef blaze::SparseColumn<MatrixType>                  ColumnType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st column of matrix A
   ColumnType col31 = column( A, 31UL );

   for( ColumnType::Iterator it=col31.begin(); it!=col31.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( ColumnType::ConstIterator it=col31.begin(); it!=col31.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_column_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse column can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedMatrix<double,blaze::columnMajor>  MatrixType;
   MatrixType A( 100UL, 10UL );  // Non-initialized 100x10 matrix

   typedef blaze::SparseColumn<MatrixType>  ColumnType;
   ColumnType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   // The subscript operator provides access to all possible elements of the sparse column,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse column, the element is inserted into the column.
   col0[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the column it is inserted into the column, if it is already contained
   // in the column its value is modified.
   col0.set( 45UL, -1.2 );

   // An alternative for inserting elements into the column is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the column.
   col0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse column is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the column and that the column's capacity is
   // large enough to hold the new element.
   col0.reserve( 10UL );
   col0.append( 51UL, -2.1 );
   \endcode

// \n \section sparse_column_common_operations Common Operations
//
// The current number of column elements can be obtained via the \c size() function, the current
// capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since columns are references to specific columns of a matrix,
// several operations are not possible on views, such as resizing and swapping:

   \code
   typedef blaze::CompressedMatrix<int,blaze::columnMajor>  MatrixType;
   typedef blaze::SparseColumn<MatrixType>                  ColumnType;

   MatrixType A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd column of matrix A
   ColumnType col2 = column( A, 2UL );

   col2.size();          // Returns the number of elements in the column
   col2.capacity();      // Returns the capacity of the column
   col2.nonZeros();      // Returns the number of non-zero elements contained in the column

   col2.resize( 84UL );  // Compilation error: Cannot resize a single column of a matrix

   ColumnType col3 = column( A, 3UL );
   swap( col2, col3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section sparse_column_arithmetic_operations Arithmetic Operations
//
// The following example gives an impression of the use of SparseColumn within arithmetic
// operations. All operations (addition, subtraction, multiplication, scaling, ...) can be
// performed on all possible combinations of dense and sparse vectors with fitting element
// types:

   \code
   blaze::CompressedVector<double,blaze::columnVector> a( 2UL ), b;
   a[1] = 2.0;
   blaze::DynamicVector<double,blaze::columnVector> c( 2UL, 3.0 );

   typedef blaze::CompressedMatrix<double,blaze::columnMajor>  MatrixType;
   MatrixType A( 2UL, 3UL );  // Non-initialized 2x3 matrix

   typedef blaze::SparseColumn<MatrixType>  ColumnType;
   ColumnType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   col0[0] = 0.0;         // Manual initialization of the 0th column of A
   col0[1] = 0.0;
   column( A, 1UL ) = a;  // Sparse vector initialization of the 1st column of A
   column( A, 2UL ) = c;  // Dense vector initialization of the 2nd column of A

   b = col0 + a;                 // Sparse vector/sparse vector addition
   b = c + column( A, 1UL );     // Dense vector/sparse vector addition
   b = col0 * column( A, 2UL );  // Component-wise vector multiplication

   column( A, 1UL ) *= 2.0;     // In-place scaling of the 1st column
   b = column( A, 1UL ) * 2.0;  // Scaling of the 1st column
   b = 2.0 * column( A, 1UL );  // Scaling of the 1st column

   column( A, 2UL ) += a;                 // Addition assignment
   column( A, 2UL ) -= c;                 // Subtraction assignment
   column( A, 2UL ) *= column( A, 0UL );  // Multiplication assignment

   double scalar = trans( c ) * column( A, 1UL );  // Scalar/dot/inner product between two vectors

   A = column( A, 1UL ) * trans( c );  // Outer product between two vectors
   \endcode

// \n \section sparse_column_on_row_major_matrix Sparse Column on a Row-Major Matrix
//
// It is especially noteworthy that column views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// also possible to traverse a column of a row-major matrix. For instance:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseColumn<MatrixType>               ColumnType;

   MatrixType A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 1st column of a row-major matrix A
   ColumnType col1 = column( A, 1UL );

   for( ColumnType::Iterator it=col1.begin(); it!=col1.end(); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a column view on a matrix stored in a row-major fashion
// can result in a considerable performance decrease in comparison to a column view on a matrix
// with column-major storage format. This is due to the non-contiguous storage of the matrix
// elements. Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two row-major matrices
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 128UL, 128UL );
   blaze::CompressedMatrix<double,blaze::rowMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th column of the multiplication between A and B ...
   blaze::CompressedVector<double,blaze::columnVector> x = column( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // A with the 15th column of the row-major matrix B.
   blaze::CompressedVector<double,blaze::columnVector> x = A * column( B, 15UL );
   \endcode

// Although Blaze performs the resulting matrix/vector multiplication as efficiently as possible
// using a column-major storage order for matrix B would result in a more efficient evaluation.
*/
template< typename MT                               // Type of the sparse matrix
        , bool SO = IsColumnMajorMatrix<MT>::value  // Storage order
        , bool SF = IsSymmetric<MT>::value >        // Symmetry flag
class SparseColumn : public SparseVector< SparseColumn<MT,SO,SF>, false >
                   , private Column
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseColumn<MT,SO,SF>              This;           //!< Type of this SparseColumn instance.
   typedef typename ColumnTrait<MT>::Type      ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the column elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseColumn&                 CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant column value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, typename MT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseColumn( MT& matrix, size_t index );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline SparseColumn& operator=( const SparseColumn& rhs );

   template< typename VT > inline SparseColumn& operator= ( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline SparseColumn& operator= ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator+=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline SparseColumn& operator+=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator-=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline SparseColumn& operator-=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator*=( const Vector<VT,false>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseColumn >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseColumn >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t        size() const;
                              inline size_t        capacity() const;
                              inline size_t        nonZeros() const;
                              inline void          reset();
                              inline Iterator      set    ( size_t index, const ElementType& value );
                              inline Iterator      insert ( size_t index, const ElementType& value );
                              inline void          erase  ( size_t index );
                              inline Iterator      erase  ( Iterator pos );
                              inline Iterator      erase  ( Iterator first, Iterator last );
                              inline void          reserve( size_t n );
   template< typename Other > inline SparseColumn& scale  ( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Low-level utility functions*****************************************************************
   /*!\name Low-level utility functions */
   //@{
   inline void append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   template< typename VT >    inline void assign   ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The sparse matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const SparseColumn<MT2,SO2,SF2>& column );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const SparseColumn<MT2,SO2,SF2>& a, const SparseColumn<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< SparseColumn<MT2,SO2,SF2> >::Type
      derestrict( SparseColumn<MT2,SO2,SF2>& dm );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
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
/*!\brief The constructor for SparseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline SparseColumn<MT,SO,SF>::SparseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the column
   , col_   ( index  )  // The index of the column in the matrix
{
   if( matrix_.columns() <= index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Reference
   SparseColumn<MT,SO,SF>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,col_);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstReference
   SparseColumn<MT,SO,SF>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,col_);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Reference
   SparseColumn<MT,SO,SF>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstReference
   SparseColumn<MT,SO,SF>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator SparseColumn<MT,SO,SF>::begin()
{
   return matrix_.begin( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator SparseColumn<MT,SO,SF>::begin() const
{
   return matrix_.cbegin( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator SparseColumn<MT,SO,SF>::cbegin() const
{
   return matrix_.cbegin( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator SparseColumn<MT,SO,SF>::end()
{
   return matrix_.end( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator SparseColumn<MT,SO,SF>::end() const
{
   return matrix_.cend( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator SparseColumn<MT,SO,SF>::cend() const
{
   return matrix_.cend( col_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for SparseColumn.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator=( const SparseColumn& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && col_ == rhs.col_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( rhs.nonZeros() );
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for dense vectors.
//
// \param rhs Dense vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      left.reset();
      assign( left, tmp );
   }
   else {
      left.reset();
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for sparse vectors.
//
// \param rhs Sparse vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( right.nonZeros() );
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a dense vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator+=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a sparse vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator+=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a dense vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator-=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a sparse vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator-=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename MultTrait<ResultType,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a sparse column
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse column.
//
// Via this operator it is possible to scale the sparse column. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for columns on lower
// or upper unitriangular matrices. The attempt to scale such a column results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse column must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,SO,SF> >::Type&
   SparseColumn<MT,SO,SF>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse column by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse column.
//
// Via this operator it is possible to scale the sparse column. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for columns on lower
// or upper unitriangular matrices. The attempt to scale such a column results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse column must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,SO,SF> >::Type&
   SparseColumn<MT,SO,SF>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( Iterator element=begin(); element!=end(); ++element )
         element->value() *= tmp;
   }
   else {
      for( Iterator element=begin(); element!=end(); ++element )
         element->value() /= rhs;
   }

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the sparse column.
//
// \return The size of the sparse column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseColumn<MT,SO,SF>::size() const
{
   return matrix_.rows();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the sparse column.
//
// \return The capacity of the sparse column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseColumn<MT,SO,SF>::capacity() const
{
   return matrix_.capacity( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseColumn<MT,SO,SF>::nonZeros() const
{
   return matrix_.nonZeros( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void SparseColumn<MT,SO,SF>::reset()
{
   matrix_.reset( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an element of the sparse column.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse column. In case the sparse column
// already contains an element with index \a index its value is modified, else a new element
// with the given \a value is inserted.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator
   SparseColumn<MT,SO,SF>::set( size_t index, const ElementType& value )
{
   return matrix_.set( index, col_, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the sparse column.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse column access index.
//
// This function inserts a new element into the sparse column. However, duplicate elements
// are not allowed. In case the sparse column already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator
   SparseColumn<MT,SO,SF>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( index, col_, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void SparseColumn<MT,SO,SF>::erase( size_t index )
{
   matrix_.erase( index, col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse column.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator SparseColumn<MT,SO,SF>::erase( Iterator pos )
{
   return matrix_.erase( col_, pos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the sparse column.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse column.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator
   SparseColumn<MT,SO,SF>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( col_, first, last );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
void SparseColumn<MT,SO,SF>::reserve( size_t n )
{
   matrix_.reserve( col_, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the scalar value
inline SparseColumn<MT,SO,SF>& SparseColumn<MT,SO,SF>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating a new sparse column capacity.
//
// \return The new sparse column capacity.
//
// This function calculates a new column capacity based on the current capacity of the sparse
// column. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseColumn<MT,SO,SF>::extendCapacity() const
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity()+1UL );
   nonzeros = max( nonzeros, 7UL    );
   nonzeros = min( nonzeros, size() );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity(), "Invalid capacity value" );

   return nonzeros;
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator SparseColumn<MT,SO,SF>::find( size_t index )
{
   return matrix_.find( index, col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator
   SparseColumn<MT,SO,SF>::find( size_t index ) const
{
   return matrix_.find( index, col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator SparseColumn<MT,SO,SF>::lowerBound( size_t index )
{
   return matrix_.lowerBound( index, col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator
   SparseColumn<MT,SO,SF>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( index, col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::Iterator SparseColumn<MT,SO,SF>::upperBound( size_t index )
{
   return matrix_.upperBound( index, col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseColumn<MT,SO,SF>::ConstIterator
   SparseColumn<MT,SO,SF>::upperBound( size_t index ) const
{
   return matrix_.upperBound( index, col_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  LOW-LEVEL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Appending an element to the sparse column.
//
// \param index The index of the new element. The index must be smaller than the number of matrix rows.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column with elements. It appends
// a new element to the end of the sparse column without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse column
//  - the current number of non-zero elements must be smaller than the capacity of the column
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void SparseColumn<MT,SO,SF>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( index, col_, value, check );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the sparse column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address can alias with the sparse column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,SO,SF>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the sparse column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,SO,SF>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,SO,SF>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<size(); ++i )
   {
      if( matrix_.nonZeros( col_ ) == matrix_.capacity( col_ ) )
         matrix_.reserve( col_, extendCapacity() );

      matrix_.append( i, col_, (~rhs)[i], true );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,SO,SF>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      matrix_.append( element->index(), col_, element->value(), true );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,SO,SF>::addAssign( const DenseVector<VT,false>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( col_ );
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,SO,SF>::addAssign( const SparseVector<VT,false>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( col_ );
   matrix_.reserve( col_, tmp.nonZeros() );
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,SO,SF>::subAssign( const DenseVector<VT,false>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( col_ );
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,SO,SF>::subAssign( const SparseVector<VT,false>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( col_ );
   matrix_.reserve( col_, tmp.nonZeros() );
   assign( tmp );
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SparseColumn for general row-major matrices.
// \ingroup views
//
// This specialization of SparseColumn adapts the class template to the requirements of general
// row-major matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SparseColumn<MT,false,false> : public SparseVector< SparseColumn<MT,false,false>, false >
                                   , private Column
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseColumn<MT,false,false>        This;            //!< Type of this SparseColumn instance.
   typedef typename ColumnTrait<MT>::Type      ResultType;      //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;     //!< Type of the column elements.
   typedef typename MT::ReturnType             ReturnType;      //!< Return type for expression template evaluations
   typedef const SparseColumn&                 CompositeType;   //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant column value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;
   //**********************************************************************************************

   //**ColumnElement class definition**************************************************************
   /*!\brief Access proxy for a specific element of the sparse column.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnElement : private SparseElement
   {
    private:
      //*******************************************************************************************
      //! Compilation switch for the return type of the value member function.
      /*! The \a returnConst compile time constant expression represents a compilation switch for
          the return type of the value member function. In case the given matrix type \a MatrixType
          is const qualified, \a returnConst will be set to 1 and the value member function will
          return a reference to const. Otherwise \a returnConst will be set to 0 and the value
          member function will offer write access to the sparse matrix elements. */
      enum { returnConst = IsConst<MatrixType>::value };
      //*******************************************************************************************

      //**Type definitions*************************************************************************
      //! Type of the underlying sparse elements.
      typedef typename std::iterator_traits<IteratorType>::value_type  SET;

      typedef typename SET::Reference       RT;   //!< Reference type of the underlying sparse element.
      typedef typename SET::ConstReference  CRT;  //!< Reference-to-const type of the underlying sparse element.
      //*******************************************************************************************

    public:
      //**Type definitions*************************************************************************
      typedef typename SET::ValueType                    ValueType;       //!< The value type of the row element.
      typedef size_t                                     IndexType;       //!< The index type of the row element.
      typedef typename IfTrue<returnConst,CRT,RT>::Type  Reference;       //!< Reference return type
      typedef CRT                                        ConstReference;  //!< Reference-to-const return type.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnElement class.
      //
      // \param pos Iterator to the current position within the sparse column.
      // \param row The row index.
      */
      inline ColumnElement( IteratorType pos, size_t row )
         : pos_( pos )  // Iterator to the current position within the sparse column
         , row_( row )  // Index of the according row
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse column element.
      //
      // \param value The new value of the sparse column element.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the addition.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the subtraction.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the multiplication.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the division.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const ColumnElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse column element.
      //
      // \return The current value of the sparse column element.
      */
      inline Reference value() const {
         return pos_->value();
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline IndexType index() const {
         return row_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse column.
      size_t       row_;  //!< Index of the according row.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**ColumnIterator class definition*************************************************************
   /*!\brief Iterator over the elements of the sparse column.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::forward_iterator_tag               IteratorCategory;  //!< The iterator category.
      typedef ColumnElement<MatrixType,IteratorType>  ValueType;         //!< Type of the underlying elements.
      typedef ValueType                               PointerType;       //!< Pointer return type.
      typedef ValueType                               ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                               DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the ColumnIterator class.
      */
      inline ColumnIterator()
         : matrix_( NULL )  // The sparse matrix containing the column.
         , row_   ( 0UL  )  // The current row index.
         , column_( 0UL  )  // The current column index.
         , pos_   ()        // Iterator to the current sparse element.
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( &matrix )  // The sparse matrix containing the column.
         , row_   ( row     )  // The current row index.
         , column_( column  )  // The current column index.
         , pos_   ()           // Iterator to the current sparse element.
      {
         for( ; row_<matrix_->rows(); ++row_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( row_ ) ) break;
         }
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      // \param pos Initial position of the iterator
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column, IteratorType pos )
         : matrix_( &matrix )  // The sparse matrix containing the column.
         , row_   ( row     )  // The current row index.
         , column_( column  )  // The current column index.
         , pos_   ( pos     )  // Iterator to the current sparse element.
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row, column ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ColumnIterator instances.
      //
      // \param it The column iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline ColumnIterator( const ColumnIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the column.
         , row_   ( it.row_    )  // The current row index.
         , column_( it.column_ )  // The current column index.
         , pos_   ( it.pos_    )  // Iterator to the current sparse element.
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ColumnIterator& operator++() {
         ++row_;
         for( ; row_<matrix_->rows(); ++row_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( row_ ) ) break;
         }

         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ColumnIterator operator++( int ) {
         const ColumnIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, row_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, row_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ == rhs.column_ );
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two column iterators.
      //
      // \param rhs The right-hand side column iterator.
      // \return The number of elements between the two column iterators.
      */
      inline DifferenceType operator-( const ColumnIterator& rhs ) const {
         size_t counter( 0UL );
         for( size_t i=rhs.row_; i<row_; ++i ) {
            if( matrix_->find( i, column_ ) != matrix_->end( i ) )
               ++counter;
         }
         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The sparse matrix containing the column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class ColumnIterator;
      template< typename MT2, bool SO2, bool SF2 > friend class SparseColumn;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef ColumnIterator<const MT,typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, ColumnIterator<MT,typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseColumn( MT& matrix, size_t index );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
                           inline SparseColumn& operator= ( const SparseColumn& rhs );
   template< typename VT > inline SparseColumn& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator*=( const Vector<VT,false>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseColumn >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseColumn >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t        size() const;
                              inline size_t        capacity() const;
                              inline size_t        nonZeros() const;
                              inline void          reset();
                              inline Iterator      set    ( size_t index, const ElementType& value );
                              inline Iterator      insert ( size_t index, const ElementType& value );
                              inline void          erase  ( size_t index );
                              inline Iterator      erase  ( Iterator pos );
                              inline Iterator      erase  ( Iterator first, Iterator last );
                              inline void          reserve( size_t n );
   template< typename Other > inline SparseColumn& scale  ( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Low-level utility functions*****************************************************************
   /*!\name Low-level utility functions */
   //@{
   inline void append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   template< typename VT >    inline void assign   ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign( const Vector<VT,false>& rhs );
   template< typename VT >    inline void subAssign( const Vector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The sparse matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const SparseColumn<MT2,SO2,SF2>& column );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const SparseColumn<MT2,SO2,SF2>& a, const SparseColumn<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< SparseColumn<MT2,SO2,SF2> >::Type
      derestrict( SparseColumn<MT2,SO2,SF2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE       ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for SparseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseColumn<MT,false,false>::SparseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the column
   , col_   ( index  )  // The index of the column in the matrix
{
   if( matrix_.columns() <= index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Reference
   SparseColumn<MT,false,false>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,col_);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstReference
   SparseColumn<MT,false,false>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,col_);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Reference
   SparseColumn<MT,false,false>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstReference
   SparseColumn<MT,false,false>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator SparseColumn<MT,false,false>::begin()
{
   return Iterator( matrix_, 0UL, col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::begin() const
{
   return ConstIterator( matrix_, 0UL, col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::cbegin() const
{
   return ConstIterator( matrix_, 0UL, col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator SparseColumn<MT,false,false>::end()
{
   return Iterator( matrix_, size(), col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::end() const
{
   return ConstIterator( matrix_, size(), col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::cend() const
{
   return ConstIterator( matrix_, size(), col_ );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for SparseColumn.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseColumn<MT,false,false>&
   SparseColumn<MT,false,false>::operator=( const SparseColumn& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && col_ == rhs.col_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      assign( left, tmp );
   }
   else {
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false,false>&
   SparseColumn<MT,false,false>::operator=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const typename VT::CompositeType tmp( ~rhs );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false,false>&
   SparseColumn<MT,false,false>::operator+=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false,false>&
   SparseColumn<MT,false,false>::operator-=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false,false>&
   SparseColumn<MT,false,false>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename MultTrait<ResultType,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a sparse column
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse column.
//
// Via this operator it is possible to scale the sparse column. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for columns on lower
// or upper unitriangular matrices. The attempt to scale such a column results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse column must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,false,false> >::Type&
   SparseColumn<MT,false,false>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= rhs;
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a sparse column by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse column.
//
// Via this operator it is possible to scale the sparse column. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for columns on lower
// or upper unitriangular matrices. The attempt to scale such a column results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse column must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,false,false> >::Type&
   SparseColumn<MT,false,false>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( Iterator element=begin(); element!=end(); ++element )
         element->value() *= tmp;
   }
   else {
      for( Iterator element=begin(); element!=end(); ++element )
         element->value() /= rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,false>::size() const
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column.
//
// \return The capacity of the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,false>::capacity() const
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,false>::nonZeros() const
{
   size_t counter( 0UL );
   for( ConstIterator element=begin(); element!=end(); ++element ) {
      ++counter;
   }
   return counter;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false,false>::reset()
{
   const size_t ibegin( ( IsLower<MT>::value )
                        ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                           ?( col_+1UL )
                           :( col_ ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper<MT>::value )
                        ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                           ?( col_ )
                           :( col_+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_.erase( i, col_ );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse column.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse column. In case the sparse column
// already contains an element with index \a index its value is modified, else a new element
// with the given \a value is inserted.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::set( size_t index, const ElementType& value )
{
   return Iterator( matrix_, index, col_, matrix_.set( index, col_, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse column access index.
//
// This function inserts a new element into the sparse column. However, duplicate elements
// are not allowed. In case the sparse column already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::insert( size_t index, const ElementType& value )
{
   return Iterator( matrix_, index, col_, matrix_.insert( index, col_, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false,false>::erase( size_t index )
{
   matrix_.erase( index, col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::erase( Iterator pos )
{
   const size_t row( pos.row_ );

   if( row == size() )
      return pos;

   matrix_.erase( row, pos.pos_ );
   return Iterator( matrix_, row+1UL, col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::erase( Iterator first, Iterator last )
{
   for( ; first!=last; ++first ) {
      matrix_.erase( first.row_, first.pos_ );
   }
   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT >  // Type of the sparse matrix
void SparseColumn<MT,false,false>::reserve( size_t n )
{
   UNUSED_PARAMETER( n );

   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the scalar value
inline SparseColumn<MT,false,false>& SparseColumn<MT,false,false>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::find( size_t index )
{
   const typename MT::Iterator pos( matrix_.find( index, col_ ) );

   if( pos != matrix_.end( index ) )
      return Iterator( matrix_, index, col_, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::find( size_t index ) const
{
   const typename MT::ConstIterator pos( matrix_.find( index, col_ ) );

   if( pos != matrix_.end( index ) )
      return ConstIterator( matrix_, index, col_, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::lowerBound( size_t index )
{
   for( size_t i=index; i<size(); ++i )
   {
      const typename MT::Iterator pos( matrix_.find( i, col_ ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, i, col_, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::lowerBound( size_t index ) const
{
   for( size_t i=index; i<size(); ++i )
   {
      const typename MT::ConstIterator pos( matrix_.find( i, col_ ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, i, col_, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::Iterator
   SparseColumn<MT,false,false>::upperBound( size_t index )
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const typename MT::Iterator pos( matrix_.find( i, col_ ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, i, col_, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,false>::ConstIterator
   SparseColumn<MT,false,false>::upperBound( size_t index ) const
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const typename MT::ConstIterator pos( matrix_.find( i, col_ ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, i, col_, pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOW-LEVEL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse column.
//
// \param index The index of the new element. The index must be smaller than the number of matrix rows.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column with elements. It appends
// a new element to the end of the sparse column without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse column
//  - the current number of non-zero elements must be smaller than the capacity of the column
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false,false>::append( size_t index, const ElementType& value, bool check )
{
   if( !check || !isDefault( value ) )
      matrix_.insert( index, col_, value );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address can alias with the sparse column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,false,false>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column, \a false if not.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,false,false>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,false,false>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      matrix_(i,col_) = (~rhs)[i];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,false,false>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      for( ; i<element->index(); ++i )
         matrix_.erase( i, col_ );
      matrix_(i++,col_) = element->value();
   }
   for( ; i<size(); ++i ) {
      matrix_.erase( i, col_ );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a vector.
//
// \param rhs The right-hand side vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline void SparseColumn<MT,false,false>::addAssign( const Vector<VT,false>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a vector.
//
// \param rhs The right-hand side vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline void SparseColumn<MT,false,false>::subAssign( const Vector<VT,false>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SparseColumn for symmetric row-major matrices.
// \ingroup views
//
// This specialization of SparseColumn adapts the class template to the requirements of symmetric
// row-major matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SparseColumn<MT,false,true> : public SparseVector< SparseColumn<MT,false,true>, false >
                                  , private Column
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseColumn<MT,false,true>         This;           //!< Type of this SparseColumn instance.
   typedef typename ColumnTrait<MT>::Type      ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the column elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseColumn&                 CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant column value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, typename MT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseColumn( MT& matrix, size_t index );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline SparseColumn& operator=( const SparseColumn& rhs );

   template< typename VT > inline SparseColumn& operator= ( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline SparseColumn& operator= ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator+=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline SparseColumn& operator+=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator-=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline SparseColumn& operator-=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator*=( const Vector<VT,false>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseColumn >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseColumn >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t        size() const;
                              inline size_t        capacity() const;
                              inline size_t        nonZeros() const;
                              inline void          reset();
                              inline Iterator      set    ( size_t index, const ElementType& value );
                              inline Iterator      insert ( size_t index, const ElementType& value );
                              inline void          erase  ( size_t index );
                              inline Iterator      erase  ( Iterator pos );
                              inline Iterator      erase  ( Iterator first, Iterator last );
                              inline void          reserve( size_t n );
   template< typename Other > inline SparseColumn& scale  ( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Low-level utility functions*****************************************************************
   /*!\name Low-level utility functions */
   //@{
   inline void append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   template< typename VT >    inline void assign   ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The sparse matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const SparseColumn<MT2,SO2,SF2>& column );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const SparseColumn<MT2,SO2,SF2>& a, const SparseColumn<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const SparseColumn<MT2,SO2,SF2>& lhs, const Vector<VT,false>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< SparseColumn<MT2,SO2,SF2> >::Type
      derestrict( SparseColumn<MT2,SO2,SF2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for SparseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseColumn<MT,false,true>::SparseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the column
   , col_   ( index  )  // The index of the column in the matrix
{
   if( matrix_.columns() <= index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Reference
   SparseColumn<MT,false,true>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(col_,index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstReference
   SparseColumn<MT,false,true>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(col_,index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Reference
   SparseColumn<MT,false,true>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstReference
   SparseColumn<MT,false,true>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator SparseColumn<MT,false,true>::begin()
{
   return matrix_.begin( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::begin() const
{
   return matrix_.cbegin( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::cbegin() const
{
   return matrix_.cbegin( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator SparseColumn<MT,false,true>::end()
{
   return matrix_.end( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::end() const
{
   return matrix_.cend( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::cend() const
{
   return matrix_.cend( col_ );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for SparseColumn.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseColumn<MT,false,true>& SparseColumn<MT,false,true>::operator=( const SparseColumn& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && col_ == rhs.col_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( rhs.nonZeros() );
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for dense vectors.
//
// \param rhs Dense vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      left.reset();
      assign( left, tmp );
   }
   else {
      left.reset();
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for sparse vectors.
//
// \param rhs Sparse vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right);
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( right.nonZeros() );
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a dense vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator+=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a sparse vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator+=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a dense vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator-=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a sparse vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator-=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false,true>&
   SparseColumn<MT,false,true>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename MultTrait<ResultType,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, col_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a sparse column
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse column.
//
// Via this operator it is possible to scale the sparse column. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for columns on lower
// or upper unitriangular matrices. The attempt to scale such a column results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse column must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,false,true> >::Type&
   SparseColumn<MT,false,true>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= rhs;
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a sparse column by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse column.
//
// Via this operator it is possible to scale the sparse column. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for columns on lower
// or upper unitriangular matrices. The attempt to scale such a column results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse column must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,false,true> >::Type&
   SparseColumn<MT,false,true>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( Iterator element=begin(); element!=end(); ++element )
         element->value() *= tmp;
   }
   else {
      for( Iterator element=begin(); element!=end(); ++element )
         element->value() /= rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the sparse column.
//
// \return The size of the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,true>::size() const
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column.
//
// \return The capacity of the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,true>::capacity() const
{
   return matrix_.capacity( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,true>::nonZeros() const
{
   return matrix_.nonZeros( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false,true>::reset()
{
   matrix_.reset( col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse column.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse column. In case the sparse column
// already contains an element with index \a index its value is modified, else a new element
// with the given \a value is inserted.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::set( size_t index, const ElementType& value )
{
   return matrix_.set( col_, index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse column access index.
//
// This function inserts a new element into the sparse column. However, duplicate elements
// are not allowed. In case the sparse column already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( col_, index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false,true>::erase( size_t index )
{
   matrix_.erase( col_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::erase( Iterator pos )
{
   return matrix_.erase( col_, pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse column.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( col_, first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT >  // Type of the sparse matrix
void SparseColumn<MT,false,true>::reserve( size_t n )
{
   matrix_.reserve( col_, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the scalar value
inline SparseColumn<MT,false,true>& SparseColumn<MT,false,true>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new sparse column capacity.
//
// \return The new sparse column capacity.
//
// This function calculates a new column capacity based on the current capacity of the sparse
// column. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseColumn<MT,false,true>::extendCapacity() const
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity()+1UL );
   nonzeros = max( nonzeros, 7UL    );
   nonzeros = min( nonzeros, size() );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity(), "Invalid capacity value" );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::find( size_t index )
{
   return matrix_.find( col_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::find( size_t index ) const
{
   return matrix_.find( col_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::lowerBound( size_t index )
{
   return matrix_.lowerBound( col_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( col_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::Iterator
   SparseColumn<MT,false,true>::upperBound( size_t index )
{
   return matrix_.upperBound( col_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false,true>::ConstIterator
   SparseColumn<MT,false,true>::upperBound( size_t index ) const
{
   return matrix_.upperBound( col_, index );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOW-LEVEL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse column.
//
// \param index The index of the new element. The index must be smaller than the number of matrix rows.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column with elements. It appends
// a new element to the end of the sparse column without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse column
//  - the current number of non-zero elements must be smaller than the capacity of the column
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false,true>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( col_, index, value, check );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address can alias with the sparse column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,false,true>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,false,true>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,false,true>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<size(); ++i )
   {
      if( matrix_.nonZeros( col_ ) == matrix_.capacity( col_ ) )
         matrix_.reserve( col_, extendCapacity() );

      matrix_.append( col_, i, (~rhs)[i], true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,false,true>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      matrix_.append( col_, element->index(), element->value(), true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,false,true>::addAssign( const DenseVector<VT,false>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( col_ );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,false,true>::addAssign( const SparseVector<VT,false>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( col_ );
   matrix_.reserve( col_, tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,false,true>::subAssign( const DenseVector<VT,false>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( col_ );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,false,true>::subAssign( const SparseVector<VT,false>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( col_ );
   matrix_.reserve( col_, tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  SPARSECOLUMN OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseColumn operators */
//@{
template< typename MT, bool SO, bool SF >
inline void reset( SparseColumn<MT,SO,SF>& column );

template< typename MT, bool SO, bool SF >
inline void clear( SparseColumn<MT,SO,SF>& column );

template< typename MT, bool SO, bool SF >
inline bool isDefault( const SparseColumn<MT,SO,SF>& column );

template< typename MT, bool SO, bool SF >
inline bool isIntact( const SparseColumn<MT,SO,SF>& column );

template< typename MT, bool SO, bool SF >
inline bool isSame( const SparseColumn<MT,SO,SF>& a, const SparseColumn<MT,SO,SF>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given sparse column.
// \ingroup sparse_column
//
// \param column The sparse column to be resetted.
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void reset( SparseColumn<MT,SO,SF>& column )
{
   column.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given sparse column.
// \ingroup sparse_column
//
// \param column The sparse column to be cleared.
// \return void
//
// Clearing a sparse column is equivalent to resetting it via the reset() function.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void clear( SparseColumn<MT,SO,SF>& column )
{
   column.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse column is in default state.
// \ingroup sparse_column
//
// \param column The sparse column to be tested for its default state.
// \return \a true in case the given column is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse column is in default state. For instance, in case
// the column is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all column elements are 0 and \a false in case any vector element is
// not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isDefault( const SparseColumn<MT,SO,SF>& column )
{
   typedef typename SparseColumn<MT,SO,SF>::ConstIterator  ConstIterator;

   const ConstIterator end( column.end() );
   for( ConstIterator element=column.begin(); element!=end; ++element )
      if( !isDefault( element->value() ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given sparse column are intact.
// \ingroup sparse_column
//
// \param column The sparse column to be tested.
// \return \a true in case the given column's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the sparse column are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::DynamicMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isIntact( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isIntact( const SparseColumn<MT,SO,SF>& column )
{
   return ( column.col_ <= column.matrix_.columns() &&
            isIntact( column.matrix_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given sparse columns represent the same observable state.
// \ingroup sparse_column
//
// \param a The first sparse column to be tested for its state.
// \param b The second sparse column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given sparse columns refer to exactly the
// same range of the same sparse matrix. In case both columns represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isSame( const SparseColumn<MT,SO,SF>& a, const SparseColumn<MT,SO,SF>& b )
{
   return ( isSame( a.matrix_, b.matrix_ ) && ( a.col_ == b.col_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a sparse column.
// \ingroup sparse_column
//
// \param lhs The target left-hand side sparse column.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool tryAssign( const SparseColumn<MT,SO,SF>& lhs, const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.matrix_, ~rhs, index, lhs.col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a sparse column.
// \ingroup sparse_column
//
// \param lhs The target left-hand side sparse column.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool tryAddAssign( const SparseColumn<MT,SO,SF>& lhs, const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.matrix_, ~rhs, index, lhs.col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a sparse column.
// \ingroup sparse_column
//
// \param lhs The target left-hand side sparse column.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool trySubAssign( const SparseColumn<MT,SO,SF>& lhs, const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.matrix_, ~rhs, index, lhs.col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a sparse column.
// \ingroup sparse_column
//
// \param lhs The target left-hand side sparse column.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool tryMultAssign( const SparseColumn<MT,SO,SF>& lhs, const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.matrix_, ~rhs, index, lhs.col_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given sparse column.
// \ingroup sparse_column
//
// \param column The sparse column to be derestricted.
// \return Sparse column without access restrictions.
//
// This function removes all restrictions on the data access to the given sparse column. It returns
// a column object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DerestrictTrait< SparseColumn<MT,SO,SF> >::Type
   derestrict( SparseColumn<MT,SO,SF>& column )
{
   typedef typename DerestrictTrait< SparseColumn<MT,SO,SF> >::Type  ReturnType;
   return ReturnType( derestrict( column.matrix_ ), column.col_ );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF >
struct IsRestricted< SparseColumn<MT,SO,SF> > : public IsTrue< IsRestricted<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DERESTRICTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF >
struct DerestrictTrait< SparseColumn<MT,SO,SF> >
{
   typedef SparseColumn< typename RemoveReference< typename DerestrictTrait<MT>::Type >::Type >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, typename T >
struct AddTrait< SparseColumn<MT,SO,SF>, T >
{
   typedef typename AddTrait< typename ColumnTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct AddTrait< T, SparseColumn<MT,SO,SF> >
{
   typedef typename AddTrait< T, typename ColumnTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool SO, bool SF, typename T >
struct SubTrait< SparseColumn<MT,SO,SF>, T >
{
   typedef typename SubTrait< typename ColumnTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct SubTrait< T, SparseColumn<MT,SO,SF> >
{
   typedef typename SubTrait< T, typename ColumnTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool SO, bool SF, typename T >
struct MultTrait< SparseColumn<MT,SO,SF>, T >
{
   typedef typename MultTrait< typename ColumnTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct MultTrait< T, SparseColumn<MT,SO,SF> >
{
   typedef typename MultTrait< T, typename ColumnTrait<MT>::Type >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CROSSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, typename T >
struct CrossTrait< SparseColumn<MT,SO,SF>, T >
{
   typedef typename CrossTrait< typename ColumnTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct CrossTrait< T, SparseColumn<MT,SO,SF> >
{
   typedef typename CrossTrait< T, typename ColumnTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool SO, bool SF, typename T >
struct DivTrait< SparseColumn<MT,SO,SF>, T >
{
   typedef typename DivTrait< typename ColumnTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct DivTrait< T, SparseColumn<MT,SO,SF> >
{
   typedef typename DivTrait< T, typename ColumnTrait<MT>::Type >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF >
struct SubvectorTrait< SparseColumn<MT,SO,SF> >
{
   typedef typename SubvectorTrait< typename SparseColumn<MT,SO,SF>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
