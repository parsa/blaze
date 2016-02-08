//=================================================================================================
/*!
//  \file blaze/math/views/SparseRow.h
//  \brief Header file for the SparseRow class template
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

#ifndef _BLAZE_MATH_VIEWS_SPARSEROW_H_
#define _BLAZE_MATH_VIEWS_SPARSEROW_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/Row.h>
#include <blaze/math/Functions.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
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
/*!\defgroup sparse_row SparseRow
// \ingroup views
*/
/*!\brief Reference to a specific row of a sparse matrix.
// \ingroup sparse_row
//
// The SparseRow template represents a reference to a specific row of a sparse matrix primitive.
// The type of the sparse matrix is specified via the first template parameter:

   \code
   template< typename MT, bool SO >
   class SparseRow;
   \endcode

//  - MT: specifies the type of the sparse matrix primitive. SparseRow can be used with every
//        sparse matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the sparse matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//  - SF: specifies whether the given matrix is a symmetric matrix or not. Also this parameter
//        doesn't have to be explicitly defined, but is automatically derived from the first
//        template parameter.
//
//
// \n \section sparse_row_setup Setup of Sparse Rows
//
// A reference to a sparse row can be created very conveniently via the \c row() function. This
// reference can be treated as any other row vector, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The reference can also be used on
// both sides of an assignment: The row can either be used as an alias to grant write access to a
// specific row of a matrix primitive on the left-hand side of an assignment or to grant read-access
// to a specific row of a matrix primitive or expression on the right-hand side of an assignment.
// The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>   SparseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   SparseMatrixType A, B;
   // ... Resizing and initialization

   // Setting the 2nd row of matrix A to x
   blaze::SparseRow<SparseMatrixType> row2 = row( A, 2UL );
   row2 = x;

   // Setting the 3rd row of matrix B to y
   row( B, 3UL ) = y;

   // Setting x to the 1st row of matrix B
   x = row( B, 1UL );

   // Setting y to the 4th row of the result of the matrix multiplication
   y = row( A * B, 4UL );
   \endcode

// \n \section sparse_row_element_access Element access
//
// A sparse row can be used like any other row vector. For instance, the elements of the sparse
// row can be directly accessed with the subscript operator.

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the 4th row of matrix A
   blaze::SparseRow<MatrixType> row4 = row( A, 4UL );

   // Setting the 1st element of the sparse row, which corresponds
   // to the 1st element in the 4th row of matrix A
   row4[1] = 2.0;
   \endcode

// The numbering of the row elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of columns of the referenced matrix. Alternatively, the elements of
// a row can be traversed via iterators. Just as with vectors, in case of non-const rows,
// \c begin() and \c end() return an Iterator, which allows a manipulation of the non-zero
// values, in case of constant rows a ConstIterator is returned:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseRow<MatrixType>                  RowType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of matrix A
   RowType row31 = row( A, 31UL );

   for( RowType::Iterator it=row31.begin(); it!=row31.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( RowType::ConstIterator it=row31.begin(); it!=row31.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_row_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse row can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A( 10UL, 100UL );  // Non-initialized 10x100 matrix

   typedef blaze::SparseRow<MatrixType>  RowType;
   RowType row0( row( A, 0UL ) );  // Reference to the 0th row of A

   // The subscript operator provides access to all possible elements of the sparse row,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse row, the element is inserted into the row.
   row0[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the row it is inserted into the row, if it is already contained in
   // the row its value is modified.
   row0.set( 45UL, -1.2 );

   // An alternative for inserting elements into the row is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the row.
   row0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse row is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the row and that the row's capacity is large
   // enough to hold the new element.
   row0.reserve( 10UL );
   row0.append( 51UL, -2.1 );
   \endcode

// \n \section sparse_row_common_operations Common Operations
//
// The current number of row elements can be obtained via the \c size() function, the current
// capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since rows are references to specific rows of a matrix,
// several operations are not possible on views, such as resizing and swapping:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseRow<MatrixType>                  RowType;

   MatrixType A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd row of matrix A
   RowType row2 = row( A, 2UL );

   row2.size();          // Returns the number of elements in the row
   row2.capacity();      // Returns the capacity of the row
   row2.nonZeros();      // Returns the number of non-zero elements contained in the row

   row2.resize( 84UL );  // Compilation error: Cannot resize a single row of a matrix

   RowType row3 = row( A, 3UL );
   swap( row2, row3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section sparse_row_arithmetic_operations Arithmetic Operations
//
// The following example gives an impression of the use of SparseRow within arithmetic operations.
// All operations (addition, subtraction, multiplication, scaling, ...) can be performed on all
// possible combinations of dense and sparse vectors with fitting element types:

   \code
   blaze::CompressedVector<double,blaze::rowVector> a( 2UL ), b;
   a[1] = 2.0;
   blaze::DynamicVector<double,blaze::rowVector> c( 2UL, 3.0 );

   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A( 3UL, 2UL );  // Non-initialized 3x2 matrix

   typedef blaze::SparseRow<MatrixType>  RowType;
   RowType row0( row( A, 0UL ) );  // Reference to the 0th row of A

   row0[0] = 0.0;        // Manual initialization of the 0th row of A
   row0[1] = 0.0;
   row( A, 1UL ) = a;    // Sparse vector initialization of the 1st row of A
   row( A, 2UL ) = c;    // Dense vector initialization of the 2nd row of A

   b = row0 + a;              // Sparse vector/sparse vector addition
   b = c + row( A, 1UL );     // Dense vector/sparse vector addition
   b = row0 * row( A, 2UL );  // Component-wise vector multiplication

   row( A, 1UL ) *= 2.0;     // In-place scaling of the 1st row
   b = row( A, 1UL ) * 2.0;  // Scaling of the 1st row
   b = 2.0 * row( A, 1UL );  // Scaling of the 1st row

   row( A, 2UL ) += a;              // Addition assignment
   row( A, 2UL ) -= c;              // Subtraction assignment
   row( A, 2UL ) *= row( A, 0UL );  // Multiplication assignment

   double scalar = row( A, 1UL ) * trans( c );  // Scalar/dot/inner product between two vectors

   A = trans( c ) * row( A, 1UL );  // Outer product between two vectors
   \endcode

// \n \section sparse_row_on_column_major_matrix Sparse Row on a Column-Major Matrix
//
// It is especially noteworthy that row views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// also possible to traverse a row of a column-major matrix. For instance:

   \code
   typedef blaze::CompressedMatrix<int,blaze::columnMajor>  MatrixType;
   typedef blaze::SparseRow<MatrixType>                     RowType;

   MatrixType A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 1st row of a column-major matrix A
   RowType row1 = row( A, 1UL );

   for( RowType::Iterator it=row1.begin(); it!=row1.end(); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a row view on a matrix stored in a column-major fashion
// can result in a considerable performance decrease in comparison to a row view on a matrix
// with row-major storage format. This is due to the non-contiguous storage of the matrix elements.
// Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two column-major matrices
   blaze::CompressedMatrix<double,blaze::columnMajor> A( 128UL, 128UL );
   blaze::CompressedMatrix<double,blaze::columnMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th row of the multiplication between A and B ...
   blaze::CompressedVector<double,blaze::rowVector> x = row( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // the 15th row of the column-major matrix A with B.
   blaze::CompressedVector<double,blaze::rowVector> x = row( A, 15UL ) * B;
   \endcode

// Although Blaze performs the resulting vector/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix A would result in a more efficient evaluation.
*/
template< typename MT                            // Type of the sparse matrix
        , bool SO = IsRowMajorMatrix<MT>::value  // Storage order
        , bool SF = IsSymmetric<MT>::value >     // Symmetry flag
class SparseRow : public SparseVector< SparseRow<MT,SO,SF>, true >
                , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseRow<MT,SO,SF>                 This;           //!< Type of this SparseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseRow&                    CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
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
   explicit inline SparseRow( MT& matrix, size_t index );
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
   inline SparseRow& operator=( const SparseRow& rhs );

   template< typename VT > inline SparseRow& operator= ( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline SparseRow& operator= ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator+=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline SparseRow& operator+=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator-=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline SparseRow& operator-=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator*=( const Vector<VT,true>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseRow >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseRow >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t     size() const;
                              inline size_t     capacity() const;
                              inline size_t     nonZeros() const;
                              inline void       reset();
                              inline Iterator   set    ( size_t index, const ElementType& value );
                              inline Iterator   insert ( size_t index, const ElementType& value );
                              inline void       erase  ( size_t index );
                              inline Iterator   erase  ( Iterator pos );
                              inline Iterator   erase  ( Iterator first, Iterator last );
                              inline void       reserve( size_t n );
   template< typename Other > inline SparseRow& scale  ( const Other& scalar );
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

   template< typename VT >    inline void assign   ( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,true>& rhs );
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
   Operand      matrix_;  //!< The sparse matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const SparseRow<MT2,SO2,SF2>& row );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const SparseRow<MT2,SO2,SF2>& a, const SparseRow<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< SparseRow<MT2,SO2,SF2> >::Type
      derestrict( SparseRow<MT2,SO2,SF2>& dm );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
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
/*!\brief The constructor for SparseRow.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline SparseRow<MT,SO,SF>::SparseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the row
   , row_   ( index  )  // The index of the row in the matrix
{
   if( matrix_.rows() <= index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Reference SparseRow<MT,SO,SF>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row_,index);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstReference
   SparseRow<MT,SO,SF>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(row_,index);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Reference SparseRow<MT,SO,SF>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstReference SparseRow<MT,SO,SF>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator SparseRow<MT,SO,SF>::begin()
{
   return matrix_.begin( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::begin() const
{
   return matrix_.cbegin( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::cbegin() const
{
   return matrix_.cbegin( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator SparseRow<MT,SO,SF>::end()
{
   return matrix_.end( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::end() const
{
   return matrix_.cend( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::cend() const
{
   return matrix_.cend( row_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for SparseRow.
//
// \param rhs Sparse row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator=( const SparseRow& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, 0UL ) ) {
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
// \return Reference to the assigned row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
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
// \return Reference to the assigned row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be added to the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator+=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be added to the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator+=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator-=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator-=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::operator*=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename MultTrait<ResultType,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication between a sparse row
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse row.
//
// Via this operator it is possible to scale the sparse row. Note however that the function is
// subject to three restrictions. First, this operator cannot be used for rows on lower or upper
// unitriangular matrices. The attempt to scale such a row results in a compilation error!
// Second, this operator can only be used for numeric data types. And third, the elements of
// the sparse row must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseRow<MT,SO,SF> >::Type&
   SparseRow<MT,SO,SF>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse row by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse row.
//
// Via this operator it is possible to scale the sparse row. Note however that the function is
// subject to three restrictions. First, this operator cannot be used for rows on lower or upper
// unitriangular matrices. The attempt to scale such a row results in a compilation error!
// Second, this operator can only be used for numeric data types. And third, the elements of
// the sparse row must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseRow<MT,SO,SF> >::Type&
   SparseRow<MT,SO,SF>::operator/=( Other rhs )
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
/*!\brief Returns the current size/dimension of the sparse row.
//
// \return The size of the sparse row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseRow<MT,SO,SF>::size() const
{
   return matrix_.columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the sparse row.
//
// \return The capacity of the sparse row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseRow<MT,SO,SF>::capacity() const
{
   return matrix_.capacity( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseRow<MT,SO,SF>::nonZeros() const
{
   return matrix_.nonZeros( row_ );
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
inline void SparseRow<MT,SO,SF>::reset()
{
   matrix_.reset( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an element of the sparse row.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse row. In case the sparse row already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator
   SparseRow<MT,SO,SF>::set( size_t index, const ElementType& value )
{
   return matrix_.set( row_, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the sparse row.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse row access index.
//
// This function inserts a new element into the sparse row. However, duplicate elements
// are not allowed. In case the sparse row already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator
   SparseRow<MT,SO,SF>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( row_, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse row.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void SparseRow<MT,SO,SF>::erase( size_t index )
{
   matrix_.erase( row_, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse row.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator SparseRow<MT,SO,SF>::erase( Iterator pos )
{
   return matrix_.erase( row_, pos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the sparse row.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse row.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator
   SparseRow<MT,SO,SF>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( row_, first, last );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the sparse row.
//
// \param n The new minimum capacity of the sparse row.
// \return void
//
// This function increases the capacity of the sparse row to at least \a n elements. The
// current values of the row elements are preserved.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
void SparseRow<MT,SO,SF>::reserve( size_t n )
{
   matrix_.reserve( row_, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the sparse row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the sparse row.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the scalar value
inline SparseRow<MT,SO,SF>& SparseRow<MT,SO,SF>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating a new sparse row capacity.
//
// \return The new sparse row capacity.
//
// This function calculates a new row capacity based on the current capacity of the sparse
// row. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t SparseRow<MT,SO,SF>::extendCapacity() const
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
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator SparseRow<MT,SO,SF>::find( size_t index )
{
   return matrix_.find( row_, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::find( size_t index ) const
{
   return matrix_.find( row_, index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator SparseRow<MT,SO,SF>::lowerBound( size_t index )
{
   return matrix_.lowerBound( row_, index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( row_, index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::Iterator SparseRow<MT,SO,SF>::upperBound( size_t index )
{
   return matrix_.upperBound( row_, index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename SparseRow<MT,SO,SF>::ConstIterator SparseRow<MT,SO,SF>::upperBound( size_t index ) const
{
   return matrix_.upperBound( row_, index );
}
//*************************************************************************************************




//=================================================================================================
//
//  LOW-LEVEL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Appending an element to the sparse row.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row with elements. It appends
// a new element to the end of the sparse row without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse row
//  - the current number of non-zero elements must be smaller than the capacity of the row
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
inline void SparseRow<MT,SO,SF>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( row_, index, value, check );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the sparse row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address can alias with the sparse row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseRow<MT,SO,SF>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the sparse row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseRow<MT,SO,SF>::isAliased( const Other* alias ) const
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
inline void SparseRow<MT,SO,SF>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t j=0UL; j<size(); ++j )
   {
      if( matrix_.nonZeros( row_ ) == matrix_.capacity( row_ ) )
         matrix_.reserve( row_, extendCapacity() );

      matrix_.append( row_, j, (~rhs)[j], true );
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
inline void SparseRow<MT,SO,SF>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      matrix_.append( row_, element->index(), element->value(), true );
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
inline void SparseRow<MT,SO,SF>::addAssign( const DenseVector<VT,true>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( row_ );
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
inline void SparseRow<MT,SO,SF>::addAssign( const SparseVector<VT,true>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( row_ );
   matrix_.reserve( row_, tmp.nonZeros() );
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
inline void SparseRow<MT,SO,SF>::subAssign( const DenseVector<VT,true>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( row_ );
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
inline void SparseRow<MT,SO,SF>::subAssign( const SparseVector<VT,true>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( row_ );
   matrix_.reserve( row_, tmp.nonZeros() );
   assign( tmp );
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SparseRow for general column-major matrices.
// \ingroup sparse_row
//
// This specialization of SparseRow adapts the class template to the requirements of general
// column-major matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SparseRow<MT,false,false> : public SparseVector< SparseRow<MT,false,false>, true >
                                , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseRow<MT,false,false>           This;            //!< Type of this SparseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;      //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;     //!< Type of the row elements.
   typedef typename MT::ReturnType             ReturnType;      //!< Return type for expression template evaluations
   typedef const SparseRow&                    CompositeType;   //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;
   //**********************************************************************************************

   //**RowElement class definition*****************************************************************
   /*!\brief Access proxy for a specific element of the sparse row.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class RowElement : private SparseElement
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
      /*!\brief Constructor for the RowElement class.
      //
      // \param pos Iterator to the current position within the sparse row.
      // \param column The column index.
      */
      inline RowElement( IteratorType pos, size_t column )
         : pos_   ( pos    )  // Iterator to the current position within the sparse row
         , column_( column )  // Index of the according column
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse row element.
      //
      // \param value The new value of the sparse row element.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the addition.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the subtraction.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the multiplication.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the division.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const RowElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse row element.
      //
      // \return The current value of the sparse row element.
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
         return column_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse row.
      size_t column_;     //!< Index of the according column.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**RowIterator class definition****************************************************************
   /*!\brief Iterator over the elements of the sparse row.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class RowIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::forward_iterator_tag            IteratorCategory;  //!< The iterator category.
      typedef RowElement<MatrixType,IteratorType>  ValueType;         //!< Type of the underlying elements.
      typedef ValueType                            PointerType;       //!< Pointer return type.
      typedef ValueType                            ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                            DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the RowIterator class.
      */
      inline RowIterator()
         : matrix_( NULL )  // The sparse matrix containing the row.
         , row_   ( 0UL  )  // The current row index.
         , column_( 0UL  )  // The current column index.
         , pos_   ()        // Iterator to the current sparse element.
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the RowIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      */
      inline RowIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( &matrix )  // The sparse matrix containing the row.
         , row_   ( row     )  // The current row index.
         , column_( column  )  // The current column index.
         , pos_   ()           // Iterator to the current sparse element.
      {
         for( ; column_<matrix_->columns(); ++column_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( column_ ) ) break;
         }
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the RowIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      // \param pos Initial position of the iterator
      */
      inline RowIterator( MatrixType& matrix, size_t row, size_t column, IteratorType pos )
         : matrix_( &matrix )  // The sparse matrix containing the row.
         , row_   ( row     )  // The current row index.
         , column_( column  )  // The current column index.
         , pos_   ( pos     )  // Iterator to the current sparse element.
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row, column ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different RowIterator instances.
      //
      // \param it The row iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline RowIterator( const RowIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the row.
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
      inline RowIterator& operator++() {
         ++column_;
         for( ; column_<matrix_->columns(); ++column_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( column_ ) ) break;
         }

         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const RowIterator operator++( int ) {
         const RowIterator tmp( *this );
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
         return ReferenceType( pos_, column_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, column_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const RowIterator<MatrixType2,IteratorType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ == rhs.column_ );
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const RowIterator<MatrixType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two row iterators.
      //
      // \param rhs The right-hand side row iterator.
      // \return The number of elements between the two row iterators.
      */
      inline DifferenceType operator-( const RowIterator& rhs ) const {
         size_t counter( 0UL );
         for( size_t j=rhs.column_; j<column_; ++j ) {
            if( matrix_->find( row_, j ) != matrix_->end( j ) )
               ++counter;
         }
         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The sparse matrix containing the row.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class RowIterator;
      template< typename MT2, bool SO2, bool SF2 > friend class SparseRow;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef RowIterator<const MT,typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, RowIterator<MT,typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseRow( MT& matrix, size_t index );
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
                           inline SparseRow& operator= ( const SparseRow& rhs );
   template< typename VT > inline SparseRow& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator*=( const Vector<VT,true>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseRow >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseRow >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t     size() const;
                              inline size_t     capacity() const;
                              inline size_t     nonZeros() const;
                              inline void       reset();
                              inline Iterator   set    ( size_t index, const ElementType& value );
                              inline Iterator   insert ( size_t index, const ElementType& value );
                              inline void       erase  ( size_t index );
                              inline Iterator   erase  ( Iterator pos );
                              inline Iterator   erase  ( Iterator first, Iterator last );
                              inline void       reserve( size_t n );
   template< typename Other > inline SparseRow& scale  ( const Other& scalar );
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

   template< typename VT >    inline void assign   ( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void addAssign( const Vector<VT,true>& rhs );
   template< typename VT >    inline void subAssign( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The sparse matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const SparseRow<MT2,SO2,SF2>& row );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const SparseRow<MT2,SO2,SF2>& a, const SparseRow<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< SparseRow<MT2,SO2,SF2> >::Type
      derestrict( SparseRow<MT2,SO2,SF2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE ( MT );
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
/*!\brief The constructor for SparseRow.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseRow<MT,false,false>::SparseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the row
   , row_   ( index  )  // The index of the row in the matrix
{
   if( matrix_.rows() <= index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Reference
   SparseRow<MT,false,false>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row_,index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstReference
   SparseRow<MT,false,false>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(row_,index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Reference
   SparseRow<MT,false,false>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstReference
   SparseRow<MT,false,false>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator SparseRow<MT,false,false>::begin()
{
   return Iterator( matrix_, row_, 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator SparseRow<MT,false,false>::begin() const
{
   return ConstIterator( matrix_, row_, 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator SparseRow<MT,false,false>::cbegin() const
{
   return ConstIterator( matrix_, row_, 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator SparseRow<MT,false,false>::end()
{
   return Iterator( matrix_, row_, size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator SparseRow<MT,false,false>::end() const
{
   return ConstIterator( matrix_, row_, size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator SparseRow<MT,false,false>::cend() const
{
   return ConstIterator( matrix_, row_, size() );
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
/*!\brief Copy assignment operator for SparseRow.
//
// \param rhs Sparse row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseRow<MT,false,false>& SparseRow<MT,false,false>::operator=( const SparseRow& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, 0UL ) ) {
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
// \return Reference to the assigned row.
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
inline SparseRow<MT,false,false>& SparseRow<MT,false,false>::operator=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const typename VT::CompositeType tmp( ~rhs );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side vector to be added to the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,false,false>& SparseRow<MT,false,false>::operator+=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,false,false>& SparseRow<MT,false,false>::operator-=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseRow<MT,false,false>& SparseRow<MT,false,false>::operator*=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename MultTrait<ResultType,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication between a sparse row
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse row.
//
// Via this operator it is possible to scale the sparse row. Note however that the function is
// subject to three restrictions. First, this operator cannot be used for rows on lower or upper
// unitriangular matrices. The attempt to scale such a row results in a compilation error!
// Second, this operator can only be used for numeric data types. And third, the elements of
// the sparse row must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseRow<MT,false,false> >::Type&
   SparseRow<MT,false,false>::operator*=( Other rhs )
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
/*!\brief Division assignment operator for the division of a sparse row by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse row.
//
// Via this operator it is possible to scale the sparse row. Note however that the function is
// subject to three restrictions. First, this operator cannot be used for rows on lower or upper
// unitriangular matrices. The attempt to scale such a row results in a compilation error!
// Second, this operator can only be used for numeric data types. And third, the elements of
// the sparse row must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseRow<MT,false,false> >::Type&
   SparseRow<MT,false,false>::operator/=( Other rhs )
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
/*!\brief Returns the current size/dimension of the row.
//
// \return The size of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,false>::size() const
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row.
//
// \return The capacity of the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,false>::capacity() const
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,false>::nonZeros() const
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
inline void SparseRow<MT,false,false>::reset()
{
   const size_t jbegin( ( IsUpper<MT>::value )
                        ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                           ?( row_+1UL )
                           :( row_ ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower<MT>::value )
                        ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                           ?( row_ )
                           :( row_+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_.erase( row_, j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse row access index.
//
// This function inserts a new element into the sparse row. However, duplicate elements
// are not allowed. In case the sparse row already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator
   SparseRow<MT,false,false>::insert( size_t index, const ElementType& value )
{
   return Iterator( matrix_, row_, index, matrix_.insert( row_, index, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse row.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse row. In case the sparse row already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator
   SparseRow<MT,false,false>::set( size_t index, const ElementType& value )
{
   return Iterator( matrix_, row_, index, matrix_.set( row_, index, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseRow<MT,false,false>::erase( size_t index )
{
   matrix_.erase( row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator SparseRow<MT,false,false>::erase( Iterator pos )
{
   const size_t column( pos.column_ );

   if( column == size() )
      return pos;

   matrix_.erase( column, pos.pos_ );
   return Iterator( matrix_, row_, column+1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator
   SparseRow<MT,false,false>::erase( Iterator first, Iterator last )
{
   for( ; first!=last; ++first ) {
      matrix_.erase( first.column_, first.pos_ );
   }
   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row.
//
// \param n The new minimum capacity of the sparse row.
// \return void
//
// This function increases the capacity of the sparse row to at least \a n elements. The
// current values of the row elements are preserved.
*/
template< typename MT >  // Type of the sparse matrix
void SparseRow<MT,false,false>::reserve( size_t n )
{
   UNUSED_PARAMETER( n );

   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the sparse row.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the scalar value
inline SparseRow<MT,false,false>& SparseRow<MT,false,false>::scale( const Other& scalar )
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
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator SparseRow<MT,false,false>::find( size_t index )
{
   const typename MT::Iterator pos( matrix_.find( row_, index ) );

   if( pos != matrix_.end( index ) )
      return Iterator( matrix_, row_, index, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator
   SparseRow<MT,false,false>::find( size_t index ) const
{
   const typename MT::ConstIterator pos( matrix_.find( row_, index ) );

   if( pos != matrix_.end( index ) )
      return ConstIterator( matrix_, row_, index, pos );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator
   SparseRow<MT,false,false>::lowerBound( size_t index )
{
   for( size_t i=index; i<size(); ++i )
   {
      const typename MT::Iterator pos( matrix_.find( row_, i ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, row_, i, pos );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator
   SparseRow<MT,false,false>::lowerBound( size_t index ) const
{
   for( size_t i=index; i<size(); ++i )
   {
      const typename MT::ConstIterator pos( matrix_.find( row_, i ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, row_, i, pos );
   }

   return end();
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::Iterator
   SparseRow<MT,false,false>::upperBound( size_t index )
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const typename MT::Iterator pos( matrix_.find( row_, i ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, row_, i, pos );
   }

   return end();
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,false>::ConstIterator
   SparseRow<MT,false,false>::upperBound( size_t index ) const
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const typename MT::ConstIterator pos( matrix_.find( row_, i ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, row_, i, pos );
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
/*!\brief Appending an element to the sparse row.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row with elements. It appends
// a new element to the end of the sparse row without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse row
//  - the current number of non-zero elements must be smaller than the capacity of the row
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
inline void SparseRow<MT,false,false>::append( size_t index, const ElementType& value, bool check )
{
   if( !check || !isDefault( value ) )
      matrix_.insert( row_, index, value );
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
/*!\brief Returns whether the sparse row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address can alias with the sparse row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseRow<MT,false,false>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row, \a false if not.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseRow<MT,false,false>::isAliased( const Other* alias ) const
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
inline void SparseRow<MT,false,false>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( size_t j=0UL; j<(~rhs).size(); ++j ) {
      matrix_(row_,j) = (~rhs)[j];
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
inline void SparseRow<MT,false,false>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   size_t j( 0UL );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      for( ; j<element->index(); ++j )
         matrix_.erase( row_, j );
      matrix_(row_,j++) = element->value();
   }
   for( ; j<size(); ++j ) {
      matrix_.erase( row_, j );
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
inline void SparseRow<MT,false,false>::addAssign( const Vector<VT,true>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
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
inline void SparseRow<MT,false,false>::subAssign( const Vector<VT,true>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SparseRow for symmetric column-major matrices.
// \ingroup sparse_row
//
// This specialization of SparseRow adapts the class template to the requirements of symmetric
// column-major matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SparseRow<MT,false,true> : public SparseVector< SparseRow<MT,false,true>, true >
                              , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseRow<MT,false,true>            This;           //!< Type of this SparseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseRow&                    CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
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
   explicit inline SparseRow( MT& matrix, size_t index );
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
   inline SparseRow& operator=( const SparseRow& rhs );

   template< typename VT > inline SparseRow& operator= ( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline SparseRow& operator= ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator+=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline SparseRow& operator+=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator-=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline SparseRow& operator-=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline SparseRow& operator*=( const Vector<VT,true>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseRow >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseRow >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t     size() const;
                              inline size_t     capacity() const;
                              inline size_t     nonZeros() const;
                              inline void       reset();
                              inline Iterator   set    ( size_t index, const ElementType& value );
                              inline Iterator   insert ( size_t index, const ElementType& value );
                              inline void       erase  ( size_t index );
                              inline Iterator   erase  ( Iterator pos );
                              inline Iterator   erase  ( Iterator first, Iterator last );
                              inline void       reserve( size_t n );
   template< typename Other > inline SparseRow& scale  ( const Other& scalar );
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

   template< typename VT >    inline void assign   ( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,true>& rhs );
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
   Operand      matrix_;  //!< The sparse matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const SparseRow<MT2,SO2,SF2>& row );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const SparseRow<MT2,SO2,SF2>& a, const SparseRow<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const SparseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< SparseRow<MT2,SO2,SF2> >::Type
      derestrict( SparseRow<MT2,SO2,SF2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
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
/*!\brief The constructor for SparseRow.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseRow<MT,false,true>::SparseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the row
   , row_   ( index  )  // The index of the row in the matrix
{
   if( matrix_.rows() <= index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Reference
   SparseRow<MT,false,true>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(index,row_);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstReference
   SparseRow<MT,false,true>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(index,row_);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Reference
   SparseRow<MT,false,true>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstReference
   SparseRow<MT,false,true>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator SparseRow<MT,false,true>::begin()
{
   return matrix_.begin( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator SparseRow<MT,false,true>::begin() const
{
   return matrix_.cbegin( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator SparseRow<MT,false,true>::cbegin() const
{
   return matrix_.cbegin( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator SparseRow<MT,false,true>::end()
{
   return matrix_.end( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator SparseRow<MT,false,true>::end() const
{
   return matrix_.cend( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator SparseRow<MT,false,true>::cend() const
{
   return matrix_.cend( row_ );
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
/*!\brief Copy assignment operator for SparseRow.
//
// \param rhs Sparse row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseRow<MT,false,true>& SparseRow<MT,false,true>::operator=( const SparseRow& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, 0UL ) ) {
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
// \return Reference to the assigned row.
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
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
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
// \return Reference to the assigned row.
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
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
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
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a dense vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be added to the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator+=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be added to the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator+=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator-=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator-=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseRow<MT,false,true>&
   SparseRow<MT,false,true>::operator*=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   typedef typename MultTrait<ResultType,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, 0UL ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication between a sparse row
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse row.
//
// Via this operator it is possible to scale the sparse row. Note however that the function is
// subject to three restrictions. First, this operator cannot be used for rows on lower or upper
// unitriangular matrices. The attempt to scale such a row results in a compilation error!
// Second, this operator can only be used for numeric data types. And third, the elements of
// the sparse row must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseRow<MT,false,true> >::Type&
   SparseRow<MT,false,true>::operator*=( Other rhs )
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
/*!\brief Division assignment operator for the division of a sparse row by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse row.
//
// Via this operator it is possible to scale the sparse row. Note however that the function is
// subject to three restrictions. First, this operator cannot be used for rows on lower or upper
// unitriangular matrices. The attempt to scale such a row results in a compilation error!
// Second, this operator can only be used for numeric data types. And third, the elements of
// the sparse row must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseRow<MT,false,true> >::Type&
   SparseRow<MT,false,true>::operator/=( Other rhs )
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
/*!\brief Returns the current size/dimension of the sparse row.
//
// \return The size of the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,true>::size() const
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row.
//
// \return The capacity of the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,true>::capacity() const
{
   return matrix_.capacity( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,true>::nonZeros() const
{
   return matrix_.nonZeros( row_ );
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
inline void SparseRow<MT,false,true>::reset()
{
   matrix_.reset( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse row.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse row. In case the sparse row already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::set( size_t index, const ElementType& value )
{
   return matrix_.set( index, row_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse row access index.
//
// This function inserts a new element into the sparse row. However, duplicate elements
// are not allowed. In case the sparse row already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( index, row_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseRow<MT,false,true>::erase( size_t index )
{
   matrix_.erase( index, row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::erase( Iterator pos )
{
   return matrix_.erase( row_, pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse row.
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( row_, first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row.
//
// \param n The new minimum capacity of the sparse row.
// \return void
//
// This function increases the capacity of the sparse row to at least \a n elements. The
// current values of the row elements are preserved.
*/
template< typename MT >  // Type of the sparse matrix
void SparseRow<MT,false,true>::reserve( size_t n )
{
   matrix_.reserve( row_, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the sparse row.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the scalar value
inline SparseRow<MT,false,true>& SparseRow<MT,false,true>::scale( const Other& scalar )
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
/*!\brief Calculating a new sparse row capacity.
//
// \return The new sparse row capacity.
//
// This function calculates a new row capacity based on the current capacity of the sparse
// row. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT >  // Type of the sparse matrix
inline size_t SparseRow<MT,false,true>::extendCapacity() const
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
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::find( size_t index )
{
   return matrix_.find( index, row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator
   SparseRow<MT,false,true>::find( size_t index ) const
{
   return matrix_.find( index, row_ );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::lowerBound( size_t index )
{
   return matrix_.lowerBound( index, row_ );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator
   SparseRow<MT,false,true>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( index, row_ );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::Iterator
   SparseRow<MT,false,true>::upperBound( size_t index )
{
   return matrix_.upperBound( index, row_ );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator or the
// insert() function!
*/
template< typename MT >  // Type of the sparse matrix
inline typename SparseRow<MT,false,true>::ConstIterator
   SparseRow<MT,false,true>::upperBound( size_t index ) const
{
   return matrix_.upperBound( index, row_ );
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
/*!\brief Appending an element to the sparse row.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row with elements. It appends
// a new element to the end of the sparse row without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse row
//  - the current number of non-zero elements must be smaller than the capacity of the row
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
inline void SparseRow<MT,false,true>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( index, row_, value, check );
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
/*!\brief Returns whether the sparse row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address can alias with the sparse row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseRow<MT,false,true>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the foreign expression
inline bool SparseRow<MT,false,true>::isAliased( const Other* alias ) const
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
inline void SparseRow<MT,false,true>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<size(); ++i )
   {
      if( matrix_.nonZeros( row_ ) == matrix_.capacity( row_ ) )
         matrix_.reserve( row_, extendCapacity() );

      matrix_.append( i, row_, (~rhs)[i], true );
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
inline void SparseRow<MT,false,true>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      matrix_.append( element->index(), row_, element->value(), true );
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
inline void SparseRow<MT,false,true>::addAssign( const DenseVector<VT,true>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( row_ );
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
inline void SparseRow<MT,false,true>::addAssign( const SparseVector<VT,true>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   matrix_.reset( row_ );
   matrix_.reserve( row_, tmp.nonZeros() );
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
inline void SparseRow<MT,false,true>::subAssign( const DenseVector<VT,true>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( row_ );
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
inline void SparseRow<MT,false,true>::subAssign( const SparseVector<VT,true>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   matrix_.reset( row_ );
   matrix_.reserve( row_, tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  SPARSEROW OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseRow operators */
//@{
template< typename MT, bool SO, bool SF >
inline void reset( SparseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline void clear( SparseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline bool isDefault( const SparseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline bool isIntact( const SparseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline bool isSame( const SparseRow<MT,SO,SF>& a, const SparseRow<MT,SO,SF>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given sparse row.
// \ingroup sparse_row
//
// \param row The sparse row to be resetted.
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void reset( SparseRow<MT,SO,SF>& row )
{
   row.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given sparse row.
// \ingroup sparse_row
//
// \param row The sparse row to be cleared.
// \return void
//
// Clearing a sparse row is equivalent to resetting it via the reset() function.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void clear( SparseRow<MT,SO,SF>& row )
{
   row.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse row is in default state.
// \ingroup sparse_row
//
// \param row The sparse row to be tested for its default state.
// \return \a true in case the given row is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse row is in default state. For instance, in case the
// row is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all row elements are 0 and \a false in case any vector element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isDefault( const SparseRow<MT,SO,SF>& row )
{
   typedef typename SparseRow<MT,SO,SF>::ConstIterator  ConstIterator;

   const ConstIterator end( row.end() );
   for( ConstIterator element=row.begin(); element!=end; ++element )
      if( !isDefault( element->value() ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given sparse row are intact.
// \ingroup sparse_row
//
// \param row The sparse row to be tested.
// \return \a true in case the given row's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the sparse row are intact, i.e. if its state
// is valid. In case the invariants are intact, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isIntact( const SparseRow<MT,SO,SF>& row )
{
   return ( row.row_ <= row.matrix_.rows() &&
            isIntact( row.matrix_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given sparse rows represent the same observable state.
// \ingroup sparse_row
//
// \param a The first sparse row to be tested for its state.
// \param b The second sparse row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given sparse rows refer to exactly the
// same range of the same sparse matrix. In case both rows represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isSame( const SparseRow<MT,SO,SF>& a, const SparseRow<MT,SO,SF>& b )
{
   return ( isSame( a.matrix_, b.matrix_ ) && ( a.row_ == b.row_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a sparse row.
// \ingroup sparse_row
//
// \param lhs The target left-hand side sparse row.
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
inline bool tryAssign( const SparseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a sparse row.
// \ingroup sparse_row
//
// \param lhs The target left-hand side sparse row.
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
inline bool tryAddAssign( const SparseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a sparse row.
// \ingroup sparse_row
//
// \param lhs The target left-hand side sparse row.
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
inline bool trySubAssign( const SparseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a sparse row.
// \ingroup sparse_row
//
// \param lhs The target left-hand side sparse row.
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
inline bool tryMultAssign( const SparseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given sparse row.
// \ingroup sparse_row
//
// \param row The sparse row to be derestricted.
// \return Sparse row without access restrictions.
//
// This function removes all restrictions on the data access to the given sparse row. It returns a
// row object that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DerestrictTrait< SparseRow<MT,SO,SF> >::Type
   derestrict( SparseRow<MT,SO,SF>& row )
{
   typedef typename DerestrictTrait< SparseRow<MT,SO,SF> >::Type  ReturnType;
   return ReturnType( derestrict( row.matrix_ ), row.row_ );
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
struct IsRestricted< SparseRow<MT,SO,SF> > : public IsTrue< IsRestricted<MT>::value >
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
struct DerestrictTrait< SparseRow<MT,SO,SF> >
{
   typedef SparseRow< typename RemoveReference< typename DerestrictTrait<MT>::Type >::Type >  Type;
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
struct AddTrait< SparseRow<MT,SO,SF>, T >
{
   typedef typename AddTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct AddTrait< T, SparseRow<MT,SO,SF> >
{
   typedef typename AddTrait< T, typename RowTrait<MT>::Type >::Type  Type;
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
struct SubTrait< SparseRow<MT,SO,SF>, T >
{
   typedef typename SubTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct SubTrait< T, SparseRow<MT,SO,SF> >
{
   typedef typename SubTrait< T, typename RowTrait<MT>::Type >::Type  Type;
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
struct MultTrait< SparseRow<MT,SO,SF>, T >
{
   typedef typename MultTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct MultTrait< T, SparseRow<MT,SO,SF> >
{
   typedef typename MultTrait< T, typename RowTrait<MT>::Type >::Type  Type;
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
struct CrossTrait< SparseRow<MT,SO,SF>, T >
{
   typedef typename CrossTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct CrossTrait< T, SparseRow<MT,SO,SF> >
{
   typedef typename CrossTrait< T, typename RowTrait<MT>::Type >::Type  Type;
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
struct DivTrait< SparseRow<MT,SO,SF>, T >
{
   typedef typename DivTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct DivTrait< T, SparseRow<MT,SO,SF> >
{
   typedef typename DivTrait< T, typename RowTrait<MT>::Type >::Type  Type;
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
struct SubvectorTrait< SparseRow<MT,SO,SF> >
{
   typedef typename SubvectorTrait< typename SparseRow<MT,SO,SF>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
