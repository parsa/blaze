//=================================================================================================
/*!
//  \file blaze/math/views/SparseSubmatrix.h
//  \brief Header file for the SparseSubmatrix class template
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

#ifndef _BLAZE_MATH_VIEWS_SPARSESUBMATRIX_H_
#define _BLAZE_MATH_VIEWS_SPARSESUBMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <vector>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/Submatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsHermitian.h>
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
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup sparse_submatrix Sparse Submatrix
// \ingroup views
*/
/*!\brief View on a specific submatrix of a sparse matrix.
// \ingroup sparse_submatrix
//
// The SparseSubmatrix template represents a view on a specific submatrix of a sparse matrix
// primitive. The type of the sparse matrix is specified via the first template parameter:

   \code
   template< typename MT, bool AF, bool SO >
   class SparseSubmatrix;
   \endcode

//  - MT: specifies the type of the sparse matrix primitive. SparseSubmatrix can be used with
//        every sparse matrix primitive, but does not work with any matrix expression type.
//  - AF: the alignment flag specifies whether the submatrix is aligned (\a blaze::aligned) or
//        unaligned (\a blaze::unaligned). The default value is \a blaze::unaligned.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the sparse matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//
//
// \n \section sparse_submatrix_setup Setup of Sparse Submatrices
//
// A view on a sparse submatrix can be created very conveniently via the \c submatrix() function:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;

   SparseMatrixType A;
   // ... Resizing and initialization

   // Creating a sparse submatrix of size 8x16, starting in row 0 and column 4
   blaze::SparseSubmatrix<SparseMatrixType> sm = submatrix( A, 0UL, 4UL, 8UL, 16UL );
   \endcode

// This view can be treated as any other sparse matrix, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The view can also be used on both
// sides of an assignment: The submatrix can either be used as an alias to grant write access to
// a specific submatrix of a sparse matrix primitive on the left-hand side of an assignment or
// to grant read-access to a specific submatrix of a sparse matrix primitive or expression on
// the right-hand side of an assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;
   typedef blaze::DynamicMatrix<double,blaze::columnMajor>  DenseMatrixType;

   SparseMatrixType A, B;
   DenseMatrixType C;
   // ... Resizing and initialization

   // Creating a sparse submatrix of size 8x4, starting in row 0 and column 2
   blaze::SparseSubmatrix<SparseMatrixType> sm = submatrix( A, 0UL, 2UL, 8UL, 4UL );

   // Setting the submatrix of A to a 8x4 submatrix of B
   sm = submatrix( B, 0UL, 0UL, 8UL, 4UL );

   // Copying the dense matrix C into another 8x4 submatrix of A
   submatrix( A, 8UL, 2UL, 8UL, 4UL ) = C;

   // Assigning part of the result of a matrix addition to the first submatrix
   sm = submatrix( B + C, 0UL, 0UL, 8UL, 4UL );
   \endcode

// \n \section sparse_submatrix_element_access Element access
//
// A sparse submatrix can be used like any other sparse matrix. For instance, the elements of the
// sparse submatrix can be directly accessed with the function call operator:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a 8x8 submatrix, starting from position (4,4)
   blaze::SparseSubmatrix<MatrixType> sm = submatrix( A, 4UL, 4UL, 8UL, 8UL );

   // Setting the element (0,0) of the submatrix, which corresponds to
   // the element at position (4,4) in matrix A
   sm(0,0) = 2.0;
   \endcode

// Alternatively, the elements of a submatrix can be traversed via (const) iterators. Just as
// with matrices, in case of non-const submatrices, \c begin() and \c end() return an Iterator,
// which allows a manipulation of the non-zero values, in case of constant submatrices a
// ConstIterator is returned:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>            SubmatrixType;

   MatrixType A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a specific submatrix of matrix A
   SubmatrixType sm = submatrix( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( SubmatrixType::Iterator it=sm.begin(0); it!=sm.end(0); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( SubmatrixType::ConstIterator it=sm.begin(1); it!=sm.end(1); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_submatrix_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse submatrix can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A( 256UL, 512UL );  // Non-initialized matrix of size 256x512

   typedef blaze::SparseSubmatrix<MatrixType>  SubmatrixType;
   SubmatrixType sm = submatrix( A, 10UL, 10UL, 16UL, 16UL );  // View on a 16x16 submatrix of A

   // The function call operator provides access to all possible elements of the sparse submatrix,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse submatrix, the element is inserted into the
   // submatrix.
   sm(2,4) = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the submatrix it is inserted into the submatrix, if it is already contained
   // in the submatrix its value is modified.
   sm.set( 2UL, 5UL, -1.2 );

   // An alternative for inserting elements into the submatrix is the \c insert() function. However,
   // it inserts the element only in case the element is not already contained in the submatrix.
   sm.insert( 2UL, 6UL, 3.7 );

   // Just as in the case of sparse matrices, elements can also be inserted via the \c append()
   // function. In case of submatrices, \c append() also requires that the appended element's
   // index is strictly larger than the currently largest non-zero index in the according row
   // or column of the submatrix and that the according row's or column's capacity is large enough
   // to hold the new element. Note however that due to the nature of a submatrix, which may be an
   // alias to the middle of a sparse matrix, the \c append() function does not work as efficiently
   // for a submatrix as it does for a matrix.
   sm.reserve( 2UL, 10UL );
   sm.append( 2UL, 10UL, -2.1 );
   \endcode

// \n \section sparse_submatrix_common_operations Common Operations
//
// The current size of the matrix, i.e. the number of rows or columns can be obtained via the
// \c rows() and \c columns() functions, the current total capacity via the \c capacity() function,
// and the number of non-zero elements via the \c nonZeros() function. However, since submatrices
// are views on a specific submatrix of a matrix, several operations are not possible on views,
// such as resizing and swapping:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>            SubmatrixType;

   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the a 8x12 submatrix of matrix A
   SubmatrixType sm = submatrix( A, 0UL, 0UL, 8UL, 12UL );

   sm.rows();      // Returns the number of rows of the submatrix
   sm.columns();   // Returns the number of columns of the submatrix
   sm.capacity();  // Returns the capacity of the submatrix
   sm.nonZeros();  // Returns the number of non-zero elements contained in the submatrix

   sm.resize( 10UL, 8UL );  // Compilation error: Cannot resize a submatrix of a matrix

   SubmatrixType sm2 = submatrix( A, 8UL, 0UL, 12UL, 8UL );
   swap( sm, sm2 );  // Compilation error: Swap operation not allowed
   \endcode

// \n \section sparse_submatrix_arithmetic_operations Arithmetic Operations
//
// The following example gives an impression of the use of SparseSubmatrix within arithmetic
// operations. All operations (addition, subtraction, multiplication, scaling, ...) can be
// performed on all possible combinations of dense and sparse matrices with fitting element
// types:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>     DenseMatrixType;
   SparseMatrixType S1, S2, S3;
   DenseMatrixType D1, D2;

   typedef blaze::StaticVector<double,8UL,blaze::columnVector>  DenseVectorType;
   DenseVectorType a, b;

   // ... Resizing and initialization

   typedef SparseSubmatrix<SparseMatrixType>  SubmatrixType;
   SubmatrixType sm = submatrix( S1, 0UL, 0UL, 8UL, 8UL );  // View on the 8x8 submatrix of matrix S1
                                                            // starting from row 0 and column 0

   submatrix( S1, 0UL, 8UL, 8UL, 8UL ) = S2;  // Sparse matrix initialization of the 8x8 submatrix
                                              // starting in row 0 and column 8
   sm = D1;                                   // Dense matrix initialization of the second 8x8 submatrix

   S3 = sm + S2;                                    // Sparse matrix/sparse matrix addition
   D2 = D1  - submatrix( S1, 8UL, 0UL, 8UL, 8UL );  // Dense matrix/sparse matrix subtraction
   S2 = sm * submatrix( S1, 8UL, 8UL, 8UL, 8UL );   // Sparse matrix/sparse matrix multiplication

   submatrix( S1, 8UL, 0UL, 8UL, 8UL ) *= 2.0;      // In-place scaling of a submatrix of S1
   S2 = submatrix( S1, 8UL, 8UL, 8UL, 8UL ) * 2.0;  // Scaling of the a submatrix of S1
   S2 = 2.0 * sm;                                   // Scaling of the a submatrix of S1

   submatrix( S1, 0UL, 8UL, 8UL, 8UL ) += S2;  // Addition assignment
   submatrix( S1, 8UL, 0UL, 8UL, 8UL ) -= D1;  // Subtraction assignment
   submatrix( S1, 8UL, 8UL, 8UL, 8UL ) *= sm;  // Multiplication assignment

   a = submatrix( S1, 4UL, 4UL, 8UL, 8UL ) * b;  // Sparse matrix/dense vector multiplication
   \endcode

// \n \section sparse_submatrix_aligned_submatrix Aligned Submatrices
//
// Usually submatrices can be defined anywhere within a matrix. They may start at any position and
// may have an arbitrary extension (only restricted by the extension of the underlying matrix).
// However, in contrast to matrices themselves, which are always properly aligned in memory and
// therefore can provide maximum performance, this means that submatrices in general have to be
// considered to be unaligned. This can be made explicit by the \a blaze::unaligned flag:

   \code
   using blaze::unaligned;

   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;

   SparseMatrixType A;
   // ... Resizing and initialization

   // Identical creations of an unaligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::SparseSubmatrix<SparseMatrixType>           sm1 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   blaze::SparseSubmatrix<SparseMatrixType>           sm2 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   blaze::SparseSubmatrix<SparseMatrixType,unaligned> sm3 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   blaze::SparseSubmatrix<SparseMatrixType,unaligned> sm4 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// All of these calls to the \c submatrix() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned submatrix. Whereas this may provide
// full flexibility in the creation of submatrices, this might result in performance restrictions
// (even in case the specified submatrix could be aligned). However, it is also possible to create
// aligned submatrices. Aligned submatrices are identical to unaligned submatrices in all aspects,
// except that they may pose additional alignment restrictions and therefore have  less flexibility
// during creation, but don't suffer from performance penalties and provide the same performance
// as the underlying matrix. Aligned submatrices are created by explicitly specifying the
// \a blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::SparseSubmatrix<SparseMatrixType,aligned> sv = submatrix<aligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// In contrast to dense submatrices, which pose several additional alignment restrictions based on
// the used element type, sparse submatrices at this time don't pose any additional restrictions.
// Therefore aligned and unaligned sparse submatrices are truly fully identical. Note however that
// this is not true for dense submatrices (see the DenseSubmatrix class description)!
//
// \n \section sparse_submatrix_on_sparse_submatrix Submatrix on Submatrix
//
// It is also possible to create a submatrix view on another submatrix. In this context it is
// important to remember that the type returned by the \c submatrix() function is the same type
// as the type of the given submatrix, since the view on a submatrix is just another view on the
// underlying sparse matrix:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>               SubmatrixType;

   MatrixType S1;

   // ... Resizing and initialization

   // Creating a submatrix view on the sparse matrix S1
   SubmatrixType sm1 = submatrix( S1, 4UL, 4UL, 8UL, 16UL );

   // Creating a submatrix view on the sparse submatrix sm1
   SubmatrixType sm2 = submatrix( sm1, 1UL, 1UL, 4UL, 8UL );
   \endcode

// \n \section sparse_submatrix_on_symmetric_matrices Submatrix on Symmetric Matrices
//
// Submatrices can also be created on symmetric matrices (see the SymmetricMatrix class template):

   \code
   using blaze::CompressedMatrix;
   using blaze::SymmetricMatrix;
   using blaze::SparseSubmatrix;

   typedef SymmetricMatrix< CompressedMatrix<int> >    SymmetricCompressedType;
   typedef SparseSubmatrix< SymmetricCompressedType >  SubmatrixType;

   // Setup of a 16x16 symmetric matrix
   SymmetricCompressedType A( 16UL );

   // Creating a sparse submatrix of size 8x12, starting in row 2 and column 4
   SubmatrixType sm = submatrix( A, 2UL, 4UL, 8UL, 12UL );
   \endcode

// It is important to note, however, that (compound) assignments to such submatrices have a
// special restriction: The symmetry of the underlying symmetric matrix must not be broken!
// Since the modification of element \f$ a_{ij} \f$ of a symmetric matrix also modifies the
// element \f$ a_{ji} \f$, the matrix to be assigned must be structured such that the symmetry
// of the symmetric matrix is preserved. Otherwise a \a std::invalid_argument exception is
// thrown:

   \code
   using blaze::CompressedMatrix;
   using blaze::SymmetricMatrix;

   // Setup of two default 4x4 symmetric matrices
   SymmetricMatrix< CompressedMatrix<int> > A1( 4 ), A2( 4 );

   // Setup of the 3x2 compressed matrix
   //
   //       ( 0 9 )
   //   B = ( 9 8 )
   //       ( 0 7 )
   //
   CompressedMatrix<int> B( 3UL, 2UL );
   B(0,0) = 1;
   B(0,1) = 2;
   B(1,0) = 3;
   B(1,1) = 4;
   B(2,1) = 5;
   B(2,2) = 6;

   // OK: Assigning B to a submatrix of A1 such that the symmetry can be preserved
   //
   //        ( 0 0 1 2 )
   //   A1 = ( 0 0 3 4 )
   //        ( 1 3 5 6 )
   //        ( 2 4 6 0 )
   //
   submatrix( A1, 0UL, 2UL, 3UL, 2UL ) = B;  // OK

   // Error: Assigning B to a submatrix of A2 such that the symmetry cannot be preserved!
   //   The elements marked with X cannot be assigned unambiguously!
   //
   //        ( 0 1 2 0 )
   //   A2 = ( 1 3 X 0 )
   //        ( 2 X 6 0 )
   //        ( 0 0 0 0 )
   //
   submatrix( A2, 0UL, 1UL, 3UL, 2UL ) = B;  // Assignment throws an exception!
   \endcode
*/
template< typename MT                                 // Type of the sparse matrix
        , bool AF = unaligned                         // Alignment flag
        , bool SO = IsColumnMajorMatrix<MT>::value >  // Storage order
class SparseSubmatrix : public SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >
                      , private Submatrix
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseSubmatrix<MT,AF,SO>           This;           //!< Type of this SparseSubmatrix instance.
   typedef typename SubmatrixTrait<MT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the submatrix elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseSubmatrix&              CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant submatrix value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;
   //**********************************************************************************************

   //**SubmatrixElement class definition***********************************************************
   /*!\brief Access proxy for a specific element of the sparse submatrix.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class SubmatrixElement : private SparseElement
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
      /*!\brief Constructor for the SubmatrixElement class.
      //
      // \param pos Iterator to the current position within the sparse submatrix.
      // \param offset The offset within the according row/column of the sparse matrix.
      */
      inline SubmatrixElement( IteratorType pos, size_t offset )
         : pos_   ( pos    )  // Iterator to the current position within the sparse submatrix
         , offset_( offset )  // Row offset within the according sparse matrix
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse submatrix element.
      //
      // \param v The new value of the sparse submatrix element.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse submatrix element.
      //
      // \param v The right-hand side value for the addition.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse submatrix element.
      //
      // \param v The right-hand side value for the subtraction.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse submatrix element.
      //
      // \param v The right-hand side value for the multiplication.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse submatrix element.
      //
      // \param v The right-hand side value for the division.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse submatrix element at the current iterator position.
      //
      // \return Reference to the sparse submatrix element at the current iterator position.
      */
      inline const SubmatrixElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse submatrix element.
      //
      // \return The current value of the sparse submatrix element.
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
         return pos_->index() - offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse submatrix.
      size_t offset_;     //!< Offset within the according row/column of the sparse matrix.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**SubmatrixIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse submatrix.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class SubmatrixIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::forward_iterator_tag                  IteratorCategory;  //!< The iterator category.
      typedef SubmatrixElement<MatrixType,IteratorType>  ValueType;         //!< Type of the underlying elements.
      typedef ValueType                                  PointerType;       //!< Pointer return type.
      typedef ValueType                                  ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                                  DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Default constructor**********************************************************************
      /*!\brief Default constructor for the SubmatrixIterator class.
      */
      inline SubmatrixIterator()
         : pos_   ()  // Iterator to the current sparse element
         , offset_()  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the SubmatrixIterator class.
      //
      // \param iterator Iterator to the current sparse element.
      // \param index The starting index within the according row/column of the sparse matrix.
      */
      inline SubmatrixIterator( IteratorType iterator, size_t index )
         : pos_   ( iterator )  // Iterator to the current sparse element
         , offset_( index    )  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline SubmatrixIterator( const SubmatrixIterator<MatrixType2,IteratorType2>& it )
         : pos_   ( it.base()   )  // Iterator to the current sparse element.
         , offset_( it.offset() )  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubmatrixIterator& operator++() {
         ++pos_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator++( int ) {
         const SubmatrixIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse submatrix element.
      //
      // \return Reference to the current sparse submatrix element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse submatrix element.
      //
      // \return Pointer to the current sparse submatrix element.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const SubmatrixIterator<MatrixType2,IteratorType2>& rhs ) const {
         return base() == rhs.base();
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const SubmatrixIterator<MatrixType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two submatrix iterators.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return The number of elements between the two submatrix iterators.
      */
      inline DifferenceType operator-( const SubmatrixIterator& rhs ) const {
         return pos_ - rhs.pos_;
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the submatrix iterator.
      //
      // \return The current position of the submatrix iterator.
      */
      inline IteratorType base() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Offset function**************************************************************************
      /*!\brief Access to the offset of the submatrix iterator.
      //
      // \return The offset of the submatrix iterator.
      */
      inline size_t offset() const {
         return offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;     //!< Iterator to the current sparse element.
      size_t       offset_;  //!< The offset of the according row/column of the sparse matrix.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubmatrixIterator<const MT,typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, SubmatrixIterator<MT,typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseSubmatrix( Operand matrix, size_t rindex, size_t cindex, size_t m, size_t n );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline SparseSubmatrix& operator=( const SparseSubmatrix& rhs );

   template< typename MT2, bool SO2 > inline SparseSubmatrix& operator= ( const Matrix<MT2,SO2>& rhs );
   template< typename MT2, bool SO2 > inline SparseSubmatrix& operator+=( const Matrix<MT2,SO2>& rhs );
   template< typename MT2, bool SO2 > inline SparseSubmatrix& operator-=( const Matrix<MT2,SO2>& rhs );
   template< typename MT2, bool SO2 > inline SparseSubmatrix& operator*=( const Matrix<MT2,SO2>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t           row() const;
                              inline size_t           rows() const;
                              inline size_t           column() const;
                              inline size_t           columns() const;
                              inline size_t           capacity() const;
                              inline size_t           capacity( size_t i ) const;
                              inline size_t           nonZeros() const;
                              inline size_t           nonZeros( size_t i ) const;
                              inline void             reset();
                              inline void             reset( size_t i );
                              inline Iterator         set( size_t i, size_t j, const ElementType& value );
                              inline Iterator         insert( size_t i, size_t j, const ElementType& value );
                              inline void             erase( size_t i, size_t j );
                              inline Iterator         erase( size_t i, Iterator pos );
                              inline Iterator         erase( size_t i, Iterator first, Iterator last );
                              inline void             reserve( size_t nonzeros );
                                     void             reserve( size_t i, size_t nonzeros );
                              inline void             trim();
                              inline void             trim( size_t i );
                              inline SparseSubmatrix& transpose();
                              inline SparseSubmatrix& ctranspose();
   template< typename Other > inline SparseSubmatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t i, size_t j );
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline Iterator      lowerBound( size_t i, size_t j );
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline Iterator      upperBound( size_t i, size_t j );
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Low-level utility functions*****************************************************************
   /*!\name Low-level utility functions */
   //@{
   inline void append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void finalize( size_t i );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool canSMPAssign() const;

   template< typename MT2, bool SO2 > inline void assign   ( const DenseMatrix<MT2,SO2>&    rhs );
   template< typename MT2 >           inline void assign   ( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 >           inline void assign   ( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2, bool SO2 > inline void addAssign( const DenseMatrix<MT2,SO2>&    rhs );
   template< typename MT2, bool SO2 > inline void addAssign( const SparseMatrix<MT2,SO2>&   rhs );
   template< typename MT2, bool SO2 > inline void subAssign( const DenseMatrix<MT2,SO2>&    rhs );
   template< typename MT2, bool SO2 > inline void subAssign( const SparseMatrix<MT2,SO2>&   rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The sparse matrix containing the submatrix.
   const size_t row_;     //!< The first row of the submatrix.
   const size_t column_;  //!< The first column of the submatrix.
   const size_t m_;       //!< The number of rows of the submatrix.
   const size_t n_;       //!< The number of columns of the submatrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< bool AF1, typename MT2, bool AF2, bool SO2 >
   friend const SparseSubmatrix<MT2,AF1,SO2>
      submatrix( const SparseSubmatrix<MT2,AF2,SO2>& sm, size_t row, size_t column, size_t m, size_t n );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isIntact( const SparseSubmatrix<MT2,AF2,SO2>& sm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const SparseSubmatrix<MT2,AF2,SO2>& a, const SparseMatrix<MT2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const SparseMatrix<MT2,SO2>& a, const SparseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const SparseSubmatrix<MT2,AF2,SO2>& a, const SparseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool tryAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                          size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename MT3, bool SO3 >
   friend bool tryAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Matrix<MT3,SO3>& rhs,
                          size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool tryAddAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename MT3, bool SO3 >
   friend bool tryAddAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Matrix<MT3,SO3>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool trySubAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename MT3, bool SO3 >
   friend bool trySubAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Matrix<MT3,SO3>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool tryMultAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                              size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2 >
   friend typename DerestrictTrait< SparseSubmatrix<MT2,AF2,SO2> >::Type
      derestrict( SparseSubmatrix<MT2,AF2,SO2>& sm );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
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
/*!\brief The constructor for SparseSubmatrix.
//
// \param matrix The sparse matrix containing the submatrix.
// \param rindex The index of the first row of the submatrix in the given sparse matrix.
// \param cindex The index of the first column of the submatrix in the given sparse matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given sparse matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline SparseSubmatrix<MT,AF,SO>::SparseSubmatrix( Operand matrix, size_t rindex, size_t cindex, size_t m, size_t n )
   : matrix_( matrix )  // The sparse matrix containing the submatrix
   , row_   ( rindex )  // The first row of the submatrix
   , column_( cindex )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
{
   if( ( row_ + m_ > matrix_.rows() ) || ( column_ + n_ > matrix_.columns() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the sparse submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Reference
   SparseSubmatrix<MT,AF,SO>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row_+i,column_+j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the sparse submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstReference
   SparseSubmatrix<MT,AF,SO>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row_+i,column_+j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Reference
   SparseSubmatrix<MT,AF,SO>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstReference
   SparseSubmatrix<MT,AF,SO>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );

   if( column_ == 0UL )
      return Iterator( matrix_.begin( i + row_ ), column_ );
   else
      return Iterator( matrix_.lowerBound( i + row_, column_ ), column_ );
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );

   if( column_ == 0UL )
      return ConstIterator( matrix_.cbegin( i + row_ ), column_ );
   else
      return ConstIterator( matrix_.lowerBound( i + row_, column_ ), column_ );
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );

   if( column_ == 0UL )
      return ConstIterator( matrix_.cbegin( i + row_ ), column_ );
   else
      return ConstIterator( matrix_.lowerBound( i + row_, column_ ), column_ );
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );

   if( matrix_.columns() == column_ + n_ )
      return Iterator( matrix_.end( i + row_ ), column_ );
   else
      return Iterator( matrix_.lowerBound( i + row_, column_ + n_ ), column_ );
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );

   if( matrix_.columns() == column_ + n_ )
      return ConstIterator( matrix_.cend( i + row_ ), column_ );
   else
      return ConstIterator( matrix_.lowerBound( i + row_, column_ + n_ ), column_ );
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );

   if( matrix_.columns() == column_ + n_ )
      return ConstIterator( matrix_.cend( i + row_ ), column_ );
   else
      return ConstIterator( matrix_.lowerBound( i + row_, column_ + n_ ), column_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for SparseSubmatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse submatrix is initialized as a copy of the given sparse submatrix. In case the
// current sizes of the two submatrices don't match, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator=( const SparseSubmatrix& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ && column_ == rhs.column_ ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Submatrix sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, column_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      left.reset();
      assign( left, tmp );
   }
   else {
      left.reset();
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   typedef typename MT2::CompositeType  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, row_, column_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( right );
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
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator+=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, column_ ) ) {
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
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator-=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, column_ ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator*=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE        ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, column_ ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication between a sparse submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse submatrix.
//
// Via this operator it is possible to scale the sparse submatrix. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for submatrices on lower
// or upper unitriangular matrices. The attempt to scale such a submatrix results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse row must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,SO> >::Type&
   SparseSubmatrix<MT,AF,SO>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<rows(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
         element->value() *= rhs;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse submatrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse submatrix.
//
// Via this operator it is possible to scale the sparse submatrix. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for submatrices on lower
// or upper unitriangular matrices. The attempt to scale such a submatrix results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse submatrix must either support the multiplication assignment operator for the
// given floating point data type or the division assignment operator for the given integral
// data type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,SO> >::Type&
   SparseSubmatrix<MT,AF,SO>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type     DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( size_t i=0UL; i<rows(); ++i ) {
         const Iterator last( end(i) );
         for( Iterator element=begin(i); element!=last; ++element )
            element->value() *= tmp;
      }
   }
   else {
      for( size_t i=0UL; i<rows(); ++i ) {
         const Iterator last( end(i) );
         for( Iterator element=begin(i); element!=last; ++element )
            element->value() /= rhs;
      }
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
/*!\brief Returns the index of the first row of the submatrix in the underlying sparse matrix.
//
// \return The index of the first row.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::row() const
{
   return row_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of rows of the sparse submatrix.
//
// \return The number of rows of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::rows() const
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the first column of the submatrix in the underlying sparse matrix.
//
// \return The index of the first column.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::column() const
{
   return column_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of columns of the sparse submatrix.
//
// \return The number of columns of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::columns() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the sparse submatrix.
//
// \return The capacity of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::capacity() const
{
   return nonZeros() + matrix_.capacity() - matrix_.nonZeros();
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::capacity( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return nonZeros( i ) + matrix_.capacity( row_+i ) - matrix_.nonZeros( row_+i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the sparse submatrix
//
// \return The number of non-zero elements in the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i )
      nonzeros += nonZeros( i );

   return nonzeros;
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
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,AF,SO>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return end(i) - begin(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,AF,SO>::reset()
{
   for( size_t i=row_; i<row_+m_; ++i )
   {
      const size_t jbegin( ( IsUpper<MT>::value )
                           ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                              ?( max( i+1UL, column_ ) )
                              :( max( i, column_ ) ) )
                           :( column_ ) );
      const size_t jend  ( ( IsLower<MT>::value )
                           ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                              ?( min( i, column_+n_ ) )
                              :( min( i+1UL, column_+n_ ) ) )
                           :( column_+n_ ) );

      matrix_.erase( i, matrix_.lowerBound( i, jbegin ), matrix_.lowerBound( i, jend ) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column to the default initial values.
//
// \param i The index of the row/column.
// \return void
//
// This function resets the values in the specified row/column to their default value. In case
// the storage order is set to \a rowMajor the function resets the values in row \a i, in case
// the storage order is set to \a columnMajor the function resets the values in column \a i.
// Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,AF,SO>::reset( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t index( row_ + i );

   const size_t jbegin( ( IsUpper<MT>::value )
                        ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                           ?( max( i+1UL, column_ ) )
                           :( max( i, column_ ) ) )
                        :( column_ ) );
   const size_t jend  ( ( IsLower<MT>::value )
                        ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                           ?( min( i, column_+n_ ) )
                           :( min( i+1UL, column_+n_ ) ) )
                        :( column_+n_ ) );

   matrix_.erase( index, matrix_.lowerBound( index, jbegin ), matrix_.lowerBound( index, jend ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an element of the sparse submatrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse submatrix. In case the sparse matrix
// already contains an element with row index \a i and column index \a j its value is modified,
// else a new element with the given \a value is inserted.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::set( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_.set( row_+i, column_+j, value ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the sparse submatrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid sparse submatrix access index.
//
// This function inserts a new element into the sparse submatrix. However, duplicate elements are
// not allowed. In case the sparse submatrix already contains an element with row index \a i and
// column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::insert( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_.insert( row_+i, column_+j, value ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse submatrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,AF,SO>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( row_ + i, column_ + j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse submatrix.
//
// \param i The row/column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse submatrix. In case the storage order is set
// to \a rowMajor the function erases an element from row \a i, in case the storage flag is set
// to \a columnMajor the function erases an element from column \a i.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::erase( size_t i, Iterator pos )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return Iterator( matrix_.erase( row_+i, pos.base() ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the sparse submatrix.
//
// \param i The row/column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of element from the sparse submatrix. In case the storage order
// is set to \a rowMajor the function erases a range of elements element from row \a i, in case
// the storage flag is set to \a columnMajor the function erases a range of elements from column
// \a i.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::erase( size_t i, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return Iterator( matrix_.erase( row_+i, first.base(), last.base() ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the sparse submatrix.
//
// \param nonzeros The new minimum capacity of the sparse submatrix.
// \return void
//
// This function increases the capacity of the sparse submatrix to at least \a nonzeros elements.
// The current values of the submatrix elements and the individual capacities of the submatrix
// rows are preserved.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,AF,SO>::reserve( size_t nonzeros )
{
   const size_t current( capacity() );

   if( nonzeros > current ) {
      matrix_.reserve( matrix_.capacity() + nonzeros - current );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of a specific row/column of the sparse submatrix.
//
// \param i The row/column index of the new element \f$[0..M-1]\f$ or \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified row/column.
// \return void
//
// This function increases the capacity of row/column \a i of the sparse submatrix to at least
// \a nonzeros elements, but not beyond the current number of columns/rows, respectively. The
// current values of the sparse submatrix and all other individual row/column capacities are
// preserved. In case the storage order is set to \a rowMajor, the function reserves capacity
// for row \a i and the index has to be in the range \f$[0..M-1]\f$. In case the storage order
// is set to \a columnMajor, the function reserves capacity for column \a i and the index has
// to be in the range \f$[0..N-1]\f$.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
void SparseSubmatrix<MT,AF,SO>::reserve( size_t i, size_t nonzeros )
{
   const size_t current( capacity( i ) );
   const size_t index  ( row_ + i );

   if( nonzeros > current ) {
      matrix_.reserve( index, matrix_.capacity( index ) + nonzeros - current );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all excessive capacity from all rows/columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all row/column-specific reserve()
// calls. The function removes all excessive capacity from all rows (in case of a rowMajor
// matrix) or columns (in case of a columnMajor matrix). Note that this function does not
// remove the overall capacity but only reduces the capacity per row/column.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
void SparseSubmatrix<MT,AF,SO>::trim()
{
   for( size_t i=0UL; i<rows(); ++i )
      trim( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all excessive capacity of a specific row/column of the sparse matrix.
//
// \param i The index of the row/column to be trimmed (\f$[0..M-1]\f$ or \f$[0..N-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a row/column-specific reserve() call.
// It removes all excessive capacity from the specified row (in case of a rowMajor matrix)
// or column (in case of a columnMajor matrix). The excessive capacity is assigned to the
// subsequent row/column.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
void SparseSubmatrix<MT,AF,SO>::trim( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   matrix_.trim( row_ + i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline SparseSubmatrix<MT,AF,SO>& SparseSubmatrix<MT,AF,SO>::transpose()
{
   using blaze::assign;

   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row_, column_ ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( trans( *this ) );
   reset();
   assign( left, tmp );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place conjugate transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline SparseSubmatrix<MT,AF,SO>& SparseSubmatrix<MT,AF,SO>::ctranspose()
{
   using blaze::assign;

   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row_, column_ ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( ctrans( *this ) );
   reset();
   assign( left, tmp );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the sparse submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the sparse submatrix.
//
// This function scales all elements of the submatrix by the given scalar value \a scalar. Note
// that the function cannot be used to scale a submatrix on a lower or upper unitriangular matrix.
// The attempt to scale such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline SparseSubmatrix<MT,AF,SO>& SparseSubmatrix<MT,AF,SO>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<rows(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
         element->value() *= scalar;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking whether there exists an overlap in the context of a symmetric matrix.
//
// \return \a true in case an overlap exists, \a false if not.
//
// This function checks if in the context of a symmetric matrix the submatrix has an overlap with
// its counterpart. In case an overlap exists, the function return \a true, otherwise it returns
// \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool SparseSubmatrix<MT,AF,SO>::hasOverlap() const
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric<MT>::value || IsHermitian<MT>::value, "Invalid matrix detected" );

   if( ( row_ + m_ <= column_ ) || ( column_ + n_ <= row_ ) )
      return false;
   else return true;
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific submatrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// submatrix. It specifically searches for the element with row index \a i and column index
// \a j. In case the element is found, the function returns an row/column iterator to the
// element. Otherwise an iterator just past the last non-zero element of row \a i or column
// \a j (the end() iterator) is returned. Note that the returned sparse submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or
// the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::find( size_t i, size_t j )
{
   const typename MT::Iterator pos( matrix_.find( row_ + i, column_ + j ) );

   if( pos != matrix_.end( row_ + i ) )
      return Iterator( pos, column_ );
   else
      return end( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific submatrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// submatrix. It specifically searches for the element with row index \a i and column index
// \a j. In case the element is found, the function returns an row/column iterator to the
// element. Otherwise an iterator just past the last non-zero element of row \a i or column
// \a j (the end() iterator) is returned. Note that the returned sparse submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or
// the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::find( size_t i, size_t j ) const
{
   const typename MT::ConstIterator pos( matrix_.find( row_ + i, column_ + j ) );

   if( pos != matrix_.end( row_ + i ) )
      return ConstIterator( pos, column_ );
   else
      return end( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index not less then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index not less then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::lowerBound( size_t i, size_t j )
{
   return Iterator( matrix_.lowerBound( row_ + i, column_ + j ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index not less then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index not less then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::lowerBound( size_t i, size_t j ) const
{
   return ConstIterator( matrix_.lowerBound( row_ + i, column_ + j ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index greater then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index greater then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::Iterator
   SparseSubmatrix<MT,AF,SO>::upperBound( size_t i, size_t j )
{
   return Iterator( matrix_.upperBound( row_ + i, column_ + j ), column_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index greater then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index greater then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,AF,SO>::ConstIterator
   SparseSubmatrix<MT,AF,SO>::upperBound( size_t i, size_t j ) const
{
   return ConstIterator( matrix_.upperBound( row_ + i, column_ + j ), column_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  LOW-LEVEL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Appending an element to the specified row/column of the sparse submatrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse submatrix with elements. It appends
// a new element to the end of the specified row/column without any additional memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row/column of the sparse submatrix
//  - the current number of non-zero elements in the submatrix must be smaller than the capacity
//    of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse submatrix:

   \code
   using blaze::rowMajor;

   typedef blaze::CompressedMatrix<double,rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>        SubmatrixType;

   MatrixType A( 42, 54 );
   SubmatrixType B = submatrix( A, 10, 10, 4, 3 );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );        // Finalizing row 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );        // Finalizing row 1
   B.finalize( 2 );        // Finalizing the empty row 2 to prepare row 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in row 3 with column index 0
   B.finalize( 3 );        // Finalizing row 3
   \endcode

// \note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,AF,SO>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( column_ + n_ == matrix_.columns() ) {
      matrix_.append( row_ + i, column_ + j, value, check );
   }
   else if( !check || !isDefault( value ) ) {
      matrix_.insert( row_ + i, column_ + j, value );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finalizing the element insertion of a row/column.
//
// \param i The index of the row/column to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a submatrix with elements.
// After completion of row/column \a i via the append() function, this function can be called to
// finalize row/column \a i and prepare the next row/column for insertion process via append().
//
// \note: Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,AF,SO>::finalize( size_t i )
{
   matrix_.trim( row_ + i );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the submatrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,AF,SO>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,AF,SO>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the submatrix can be used in SMP assignments.
//
// \return \a true in case the submatrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the submatrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool SparseSubmatrix<MT,AF,SO>::canSMPAssign() const
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,AF,SO>::assign( const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   reserve( 0UL, rows() * columns() );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( IsSymmetric<MT>::value || IsHermitian<MT>::value )
            set( i, j, (~rhs)(i,j) );
         else
            append( i, j, (~rhs)(i,j), true );
      }
      finalize( i );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the sparse matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,AF,SO>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   reserve( 0UL, (~rhs).nonZeros() );

   for( size_t i=0UL; i<(~rhs).rows(); ++i ) {
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element ) {
         if( IsSymmetric<MT>::value || IsHermitian<MT>::value )
            set( i, element->index(), element->value() );
         else
            append( i, element->index(), element->value(), true );
      }
      finalize( i );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the sparse matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,AF,SO>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   typedef typename MT2::ConstIterator  RhsIterator;

   // Counting the number of elements per row
   std::vector<size_t> rowLengths( m_, 0UL );
   for( size_t j=0UL; j<n_; ++j ) {
      for( RhsIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         ++rowLengths[element->index()];
   }

   // Resizing the sparse matrix
   for( size_t i=0UL; i<m_; ++i ) {
      reserve( i, rowLengths[i] );
   }

   // Appending the elements to the rows of the sparse submatrix
   for( size_t j=0UL; j<n_; ++j ) {
      for( RhsIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         if( IsSymmetric<MT>::value || IsHermitian<MT>::value )
            set( element->index(), j, element->value() );
         else
            append( element->index(), j, element->value(), true );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,AF,SO>::addAssign( const DenseMatrix<MT2,SO2>& rhs )
{
   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const AddType tmp( serial( *this + (~rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,AF,SO>::addAssign( const SparseMatrix<MT2,SO2>& rhs )
{
   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const AddType tmp( serial( *this + (~rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,AF,SO>::subAssign( const DenseMatrix<MT2,SO2>& rhs )
{
   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const SubType tmp( serial( *this - (~rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand sparse matrix
inline void SparseSubmatrix<MT,AF,SO>::subAssign( const SparseMatrix<MT2,SO2>& rhs )
{
   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const SubType tmp( serial( *this - (~rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SparseSubmatrix for column-major matrices.
// \ingroup sparse_submatrix
//
// This specialization of SparseSubmatrix adapts the class template to the requirements of
// column-major matrices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
class SparseSubmatrix<MT,AF,true> : public SparseMatrix< SparseSubmatrix<MT,AF,true>, true >
                                  , private Submatrix
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseSubmatrix<MT,AF,true>         This;           //!< Type of this SparseSubmatrix instance.
   typedef typename SubmatrixTrait<MT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the submatrix elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseSubmatrix&              CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant submatrix value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;
   //**********************************************************************************************

   //**SubmatrixElement class definition***********************************************************
   /*!\brief Access proxy for a specific element of the sparse submatrix.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class SubmatrixElement : private SparseElement
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
      /*!\brief Constructor for the SubmatrixElement class.
      //
      // \param pos Iterator to the current position within the sparse submatrix.
      // \param offset The offset within the according row/column of the sparse matrix.
      */
      inline SubmatrixElement( IteratorType pos, size_t offset )
         : pos_   ( pos    )  // Iterator to the current position within the sparse submatrix
         , offset_( offset )  // Row offset within the according sparse matrix
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse submatrix element.
      //
      // \param value The new value of the sparse submatrix element.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse submatrix element.
      //
      // \param value The right-hand side value for the addition.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse submatrix element.
      //
      // \param value The right-hand side value for the subtraction.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse submatrix element.
      //
      // \param value The right-hand side value for the multiplication.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse submatrix element.
      //
      // \param value The right-hand side value for the division.
      // \return Reference to the sparse submatrix element.
      */
      template< typename T > inline SubmatrixElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse submatrix element at the current iterator position.
      //
      // \return Reference to the sparse submatrix element at the current iterator position.
      */
      inline const SubmatrixElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse submatrix element.
      //
      // \return The current value of the sparse submatrix element.
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
         return pos_->index() - offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse submatrix.
      size_t offset_;     //!< Offset within the according row/column of the sparse matrix.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**SubmatrixIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse submatrix.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class SubmatrixIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::forward_iterator_tag                  IteratorCategory;  //!< The iterator category.
      typedef SubmatrixElement<MatrixType,IteratorType>  ValueType;         //!< Type of the underlying elements.
      typedef ValueType                                  PointerType;       //!< Pointer return type.
      typedef ValueType                                  ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                                  DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Default constructor**********************************************************************
      /*!\brief Default constructor for the SubmatrixIterator class.
      */
      inline SubmatrixIterator()
         : pos_   ()  // Iterator to the current sparse element
         , offset_()  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the SubmatrixIterator class.
      //
      // \param iterator Iterator to the current sparse element.
      // \param index The starting index within the according row/column of the sparse matrix.
      */
      inline SubmatrixIterator( IteratorType iterator, size_t index )
         : pos_   ( iterator )  // Iterator to the current sparse element
         , offset_( index    )  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline SubmatrixIterator( const SubmatrixIterator<MatrixType2,IteratorType2>& it )
         : pos_   ( it.base()   )  // Iterator to the current sparse element.
         , offset_( it.offset() )  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubmatrixIterator& operator++() {
         ++pos_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator++( int ) {
         const SubmatrixIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse submatrix element.
      //
      // \return Reference to the current sparse submatrix element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse submatrix element.
      //
      // \return Pointer to the current sparse submatrix element.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const SubmatrixIterator<MatrixType2,IteratorType2>& rhs ) const {
         return base() == rhs.base();
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const SubmatrixIterator<MatrixType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two submatrix iterators.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return The number of elements between the two submatrix iterators.
      */
      inline DifferenceType operator-( const SubmatrixIterator& rhs ) const {
         return pos_ - rhs.pos_;
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the submatrix iterator.
      //
      // \return The current position of the submatrix iterator.
      */
      inline IteratorType base() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Offset function**************************************************************************
      /*!\brief Access to the offset of the submatrix iterator.
      //
      // \return The offset of the submatrix iterator.
      */
      inline size_t offset() const {
         return offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;     //!< Iterator to the current sparse element.
      size_t       offset_;  //!< The offset of the according row/column of the sparse matrix.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubmatrixIterator<const MT,typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, SubmatrixIterator<MT,typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseSubmatrix( Operand matrix, size_t rindex, size_t cindex, size_t m, size_t n );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline SparseSubmatrix& operator=( const SparseSubmatrix& rhs );

   template< typename MT2, bool SO > inline SparseSubmatrix& operator= ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline SparseSubmatrix& operator+=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline SparseSubmatrix& operator-=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline SparseSubmatrix& operator*=( const Matrix<MT2,SO>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t           row() const;
                              inline size_t           rows() const;
                              inline size_t           column() const;
                              inline size_t           columns() const;
                              inline size_t           capacity() const;
                              inline size_t           capacity( size_t i ) const;
                              inline size_t           nonZeros() const;
                              inline size_t           nonZeros( size_t i ) const;
                              inline void             reset();
                              inline void             reset( size_t i );
                              inline Iterator         set( size_t i, size_t j, const ElementType& value );
                              inline Iterator         insert( size_t i, size_t j, const ElementType& value );
                              inline void             erase( size_t i, size_t j );
                              inline Iterator         erase( size_t i, Iterator pos );
                              inline Iterator         erase( size_t i, Iterator first, Iterator last );
                              inline void             reserve( size_t nonzeros );
                                     void             reserve( size_t i, size_t nonzeros );
                              inline void             trim();
                              inline void             trim( size_t j );
                              inline SparseSubmatrix& transpose();
                              inline SparseSubmatrix& ctranspose();
   template< typename Other > inline SparseSubmatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t i, size_t j );
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline Iterator      lowerBound( size_t i, size_t j );
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline Iterator      upperBound( size_t i, size_t j );
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Low-level utility functions*****************************************************************
   /*!\name Low-level utility functions */
   //@{
   inline void append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void finalize( size_t i );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool canSMPAssign() const;

   template< typename MT2, bool SO > inline void assign   ( const DenseMatrix<MT2,SO>&     rhs );
   template< typename MT2 >          inline void assign   ( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 >          inline void assign   ( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2, bool SO > inline void addAssign( const DenseMatrix<MT2,SO>&     rhs );
   template< typename MT2, bool SO > inline void addAssign( const SparseMatrix<MT2,SO>&    rhs );
   template< typename MT2, bool SO > inline void subAssign( const DenseMatrix<MT2,SO>&     rhs );
   template< typename MT2, bool SO > inline void subAssign( const SparseMatrix<MT2,SO>&    rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The sparse matrix containing the submatrix.
   const size_t row_;     //!< The first row of the submatrix.
   const size_t column_;  //!< The first column of the submatrix.
   const size_t m_;       //!< The number of rows of the submatrix.
   const size_t n_;       //!< The number of columns of the submatrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< bool AF1, typename MT2, bool AF2, bool SO2 >
   friend const SparseSubmatrix<MT2,AF1,SO2>
      submatrix( const SparseSubmatrix<MT2,AF2,SO2>& sm, size_t row, size_t column, size_t m, size_t n );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isIntact( const SparseSubmatrix<MT2,AF2,SO2>& sm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const SparseSubmatrix<MT2,AF2,SO2>& a, const SparseMatrix<MT2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const SparseMatrix<MT2,SO2>& a, const SparseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const SparseSubmatrix<MT2,AF2,SO2>& a, const SparseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool tryAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                          size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename MT3, bool SO3 >
   friend bool tryAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Matrix<MT3,SO3>& rhs,
                          size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool tryAddAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename MT3, bool SO3 >
   friend bool tryAddAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Matrix<MT3,SO3>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool trySubAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename MT3, bool SO3 >
   friend bool trySubAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Matrix<MT3,SO3>& rhs,
                             size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2, typename VT, bool TF >
   friend bool tryMultAssign( const SparseSubmatrix<MT2,AF2,SO2>& lhs, const Vector<VT,TF>& rhs,
                              size_t row, size_t column );

   template< typename MT2, bool AF2, bool SO2 >
   friend typename DerestrictTrait< SparseSubmatrix<MT2,AF2,SO2> >::Type
      derestrict( SparseSubmatrix<MT2,AF2,SO2>& sm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
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
/*!\brief The constructor for SparseSubmatrix.
//
// \param matrix The sparse matrix containing the submatrix.
// \param rindex The index of the first row of the submatrix in the given sparse matrix.
// \param cindex The index of the first column of the submatrix in the given sparse matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given sparse matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline SparseSubmatrix<MT,AF,true>::SparseSubmatrix( Operand matrix, size_t rindex, size_t cindex, size_t m, size_t n )
   : matrix_( matrix )  // The sparse matrix containing the submatrix
   , row_   ( rindex )  // The first row of the submatrix
   , column_( cindex )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
{
   if( ( row_ + m_ > matrix_.rows() ) || ( column_ + n_ > matrix_.columns() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
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
/*!\brief 2D-access to the sparse submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Reference
   SparseSubmatrix<MT,AF,true>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the sparse submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstReference
   SparseSubmatrix<MT,AF,true>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Reference
   SparseSubmatrix<MT,AF,true>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstReference
   SparseSubmatrix<MT,AF,true>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid sparse submatrix column access index" );

   if( row_ == 0UL )
      return Iterator( matrix_.begin( j + column_ ), row_ );
   else
      return Iterator( matrix_.lowerBound( row_, j + column_ ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid sparse submatrix column access index" );

   if( row_ == 0UL )
      return ConstIterator( matrix_.cbegin( j + column_ ), row_ );
   else
      return ConstIterator( matrix_.lowerBound( row_, j + column_ ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid sparse submatrix column access index" );

   if( row_ == 0UL )
      return ConstIterator( matrix_.cbegin( j + column_ ), row_ );
   else
      return ConstIterator( matrix_.lowerBound( row_, j + column_ ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
  SparseSubmatrix<MT,AF,true>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid sparse submatrix column access index" );

   if( matrix_.rows() == row_ + m_ )
      return Iterator( matrix_.end( j + column_ ), row_ );
   else
      return Iterator( matrix_.lowerBound( row_ + m_, j + column_ ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid sparse submatrix column access index" );

   if( matrix_.rows() == row_ + m_ )
      return ConstIterator( matrix_.cend( j + column_ ), row_ );
   else
      return ConstIterator( matrix_.lowerBound( row_ + m_, j + column_ ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid sparse submatrix column access index" );

   if( matrix_.rows() == row_ + m_ )
      return ConstIterator( matrix_.cend( j + column_ ), row_ );
   else
      return ConstIterator( matrix_.lowerBound( row_ + m_, j + column_ ), row_ );
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
/*!\brief Copy assignment operator for SparseSubmatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse submatrix is initialized as a copy of the given sparse submatrix. In case the
// current sizes of the two submatrices don't match, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator=( const SparseSubmatrix& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ && column_ == rhs.column_ ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Submatrix sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, column_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      left.reset();
      assign( left, tmp );
   }
   else {
      left.reset();
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator=( const Matrix<MT2,SO>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   typedef typename MT2::CompositeType  Right;
   Right right( ~rhs );

   if( !tryAssign( matrix_, right, row_, column_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( right );
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
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator+=( const Matrix<MT2,SO>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, column_ ) ) {
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
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator-=( const Matrix<MT2,SO>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, column_ ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator*=( const Matrix<MT2,SO>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE        ( MultType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType   );

   if( columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( matrix_, tmp, row_, column_ ) ) {
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
/*!\brief Multiplication assignment operator for the multiplication between a sparse submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse submatrix.
//
// Via this operator it is possible to scale the sparse submatrix. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for submatrices on lower
// or upper unitriangular matrices. The attempt to scale such a submatrix results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse row must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,true> >::Type&
   SparseSubmatrix<MT,AF,true>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<columns(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
         element->value() *= rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a sparse submatrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the sparse submatrix.
//
// Via this operator it is possible to scale the sparse submatrix. Note however that the function
// is subject to three restrictions. First, this operator cannot be used for submatrices on lower
// or upper unitriangular matrices. The attempt to scale such a submatrix results in a compilation
// error! Second, this operator can only be used for numeric data types. And third, the elements
// of the sparse submatrix must either support the multiplication assignment operator for the
// given floating point data type or the division assignment operator for the given integral
// data type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,true> >::Type&
   SparseSubmatrix<MT,AF,true>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type     DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( size_t i=0UL; i<columns(); ++i ) {
         const Iterator last( end(i) );
         for( Iterator element=begin(i); element!=last; ++element )
            element->value() *= tmp;
      }
   }
   else {
      for( size_t i=0UL; i<columns(); ++i ) {
         const Iterator last( end(i) );
         for( Iterator element=begin(i); element!=last; ++element )
            element->value() /= rhs;
      }
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
/*!\brief Returns the index of the first row of the submatrix in the underlying sparse matrix.
//
// \return The index of the first row.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::row() const
{
   return row_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the sparse submatrix.
//
// \return The number of rows of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::rows() const
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the first column of the submatrix in the underlying sparse matrix.
//
// \return The index of the first column.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::column() const
{
   return column_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the sparse submatrix.
//
// \return The number of columns of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::columns() const
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse submatrix.
//
// \return The capacity of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::capacity() const
{
   return nonZeros() + matrix_.capacity() - matrix_.nonZeros();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::capacity( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return nonZeros( j ) + matrix_.capacity( column_+j ) - matrix_.nonZeros( column_+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse submatrix
//
// \return The number of non-zero elements in the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<columns(); ++i )
      nonzeros += nonZeros( i );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline size_t SparseSubmatrix<MT,AF,true>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return end(j) - begin(j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline void SparseSubmatrix<MT,AF,true>::reset()
{
   for( size_t j=column_; j<column_+n_; ++j )
   {
      const size_t ibegin( ( IsLower<MT>::value )
                           ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                              ?( max( j+1UL, row_ ) )
                              :( max( j, row_ ) ) )
                           :( row_ ) );
      const size_t iend  ( ( IsUpper<MT>::value )
                           ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                              ?( min( j, row_+m_ ) )
                              :( min( j+1UL, row_+m_ ) ) )
                           :( row_+m_ ) );

      matrix_.erase( j, matrix_.lowerBound( ibegin, j ), matrix_.lowerBound( iend, j ) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline void SparseSubmatrix<MT,AF,true>::reset( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t index( column_ + j );

   const size_t ibegin( ( IsLower<MT>::value )
                        ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                           ?( max( j+1UL, row_ ) )
                           :( max( j, row_ ) ) )
                        :( row_ ) );
   const size_t iend  ( ( IsUpper<MT>::value )
                        ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                           ?( min( j, row_+m_ ) )
                           :( min( j+1UL, row_+m_ ) ) )
                        :( row_+m_ ) );

   matrix_.erase( index, matrix_.lowerBound( ibegin, index ), matrix_.lowerBound( iend, index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse submatrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse submatrix. In case the sparse matrix
// already contains an element with row index \a i and column index \a j its value is modified,
// else a new element with the given \a value is inserted.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::set( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_.set( row_+i, column_+j, value ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse submatrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid sparse submatrix access index.
//
// This function inserts a new element into the sparse submatrix. However, duplicate elements are
// not allowed. In case the sparse submatrix already contains an element with row index \a i and
// column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::insert( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_.insert( row_+i, column_+j, value ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse submatrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline void SparseSubmatrix<MT,AF,true>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( row_ + i, column_ + j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse submatrix.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from column \a j of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::erase( size_t j, Iterator pos )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return Iterator( matrix_.erase( column_+j, pos.base() ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse submatrix.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of element from column \a j of the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::erase( size_t j, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return Iterator( matrix_.erase( column_+j, first.base(), last.base() ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse submatrix.
//
// \param nonzeros The new minimum capacity of the sparse submatrix.
// \return void
//
// This function increases the capacity of the sparse submatrix to at least \a nonzeros elements.
// The current values of the submatrix elements and the individual capacities of the submatrix
// rows are preserved.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline void SparseSubmatrix<MT,AF,true>::reserve( size_t nonzeros )
{
   const size_t current( capacity() );

   if( nonzeros > current ) {
      matrix_.reserve( matrix_.capacity() + nonzeros - current );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of a specific column of the sparse submatrix.
//
// \param j The column index of the new element \f$[0..M-1]\f$ or \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified column.
// \return void
//
// This function increases the capacity of column \a i of the sparse submatrix to at least
// \a nonzeros elements, but not beyond the current number of rows. The current values of
// the sparse submatrix and all other individual row/column capacities are preserved.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
void SparseSubmatrix<MT,AF,true>::reserve( size_t j, size_t nonzeros )
{
   const size_t current( capacity( j ) );
   const size_t index  ( column_ + j );

   if( nonzeros > current ) {
      matrix_.reserve( index, matrix_.capacity( index ) + nonzeros - current );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all column-specific reserve() calls
// It removes all excessive capacity from all columns. Note that this function does not remove
// the overall capacity but only reduces the capacity per column.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
void SparseSubmatrix<MT,AF,true>::trim()
{
   for( size_t j=0UL; j<columns(); ++j )
      trim( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific column of the sparse matrix.
//
// \param j The index of the column to be trimmed (\f$[0..N-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a column-specific reserve() call. It
// removes all excessive capacity from the specified column. The excessive capacity is assigned
// to the subsequent column.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
void SparseSubmatrix<MT,AF,true>::trim( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   matrix_.trim( column_ + j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline SparseSubmatrix<MT,AF,true>& SparseSubmatrix<MT,AF,true>::transpose()
{
   using blaze::assign;

   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row_, column_ ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( trans( *this ) );
   reset();
   assign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline SparseSubmatrix<MT,AF,true>& SparseSubmatrix<MT,AF,true>::ctranspose()
{
   using blaze::assign;

   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, ctrans( *this ), row_, column_ ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( ctrans(*this) );
   reset();
   assign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the sparse submatrix.
//
// This function scales all elements of the submatrix by the given scalar value \a scalar. Note
// that the function cannot be used to scale a submatrix on a lower or upper unitriangular matrix.
// The attempt to scale such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the scalar value
inline SparseSubmatrix<MT,AF,true>& SparseSubmatrix<MT,AF,true>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<columns(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
         element->value() *= scalar;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking whether there exists an overlap in the context of a symmetric matrix.
//
// \return \a true in case an overlap exists, \a false if not.
//
// This function checks if in the context of a symmetric matrix the submatrix has an overlap with
// its counterpart. In case an overlap exists, the function return \a true, otherwise it returns
// \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline bool SparseSubmatrix<MT,AF,true>::hasOverlap() const
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric<MT>::value || IsHermitian<MT>::value, "Invalid matrix detected" );

   if( ( row_ + m_ <= column_ ) || ( column_ + n_ <= row_ ) )
      return false;
   else return true;
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
/*!\brief Searches for a specific submatrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// submatrix. It specifically searches for the element with row index \a i and column index
// \a j. In case the element is found, the function returns an row/column iterator to the
// element. Otherwise an iterator just past the last non-zero element of row \a i or column
// \a j (the end() iterator) is returned. Note that the returned sparse submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or
// the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::find( size_t i, size_t j )
{
   const typename MT::Iterator pos( matrix_.find( row_ + i, column_ + j ) );

   if( pos != matrix_.end( column_ + j ) )
      return Iterator( pos, row_ );
   else
      return end( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific submatrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// submatrix. It specifically searches for the element with row index \a i and column index
// \a j. In case the element is found, the function returns an row/column iterator to the
// element. Otherwise an iterator just past the last non-zero element of row \a i or column
// \a j (the end() iterator) is returned. Note that the returned sparse submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or
// the insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::find( size_t i, size_t j ) const
{
   const typename MT::ConstIterator pos( matrix_.find( row_ + i, column_ + j ) );

   if( pos != matrix_.end( column_ + j ) )
      return ConstIterator( pos, row_ );
   else
      return end( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index not less then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index not less then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::lowerBound( size_t i, size_t j )
{
   return Iterator( matrix_.lowerBound( row_ + i, column_ + j ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index not less then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index not less then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::lowerBound( size_t i, size_t j ) const
{
   return ConstIterator( matrix_.lowerBound( row_ + i, column_ + j ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index greater then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index greater then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::Iterator
   SparseSubmatrix<MT,AF,true>::upperBound( size_t i, size_t j )
{
   return Iterator( matrix_.upperBound( row_ + i, column_ + j ), row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case of a row-major submatrix, this function returns a row iterator to the first element
// with an index greater then the given column index. In case of a column-major submatrix, the
// function returns a column iterator to the first element with an index greater then the given
// row index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned submatrix iterator
// is subject to invalidation due to inserting operations via the function call operator or the
// insert() function!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline typename SparseSubmatrix<MT,AF,true>::ConstIterator
   SparseSubmatrix<MT,AF,true>::upperBound( size_t i, size_t j ) const
{
   return ConstIterator( matrix_.upperBound( row_ + i, column_ + j ), row_ );
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
/*!\brief Appending an element to the specified row/column of the sparse submatrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse submatrix with elements. It appends
// a new element to the end of the specified row/column without any additional memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row/column of the sparse submatrix
//  - the current number of non-zero elements in the submatrix must be smaller than the capacity
//    of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse submatrix:

   \code
   using blaze::rowMajor;

   typedef blaze::CompressedMatrix<double,rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>        SubmatrixType;

   MatrixType A( 42, 54 );
   SubmatrixType B = submatrix( A, 10, 10, 4, 3 );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );        // Finalizing row 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );        // Finalizing row 1
   B.finalize( 2 );        // Finalizing the empty row 2 to prepare row 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in row 3 with column index 0
   B.finalize( 3 );        // Finalizing row 3
   \endcode

// \note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline void SparseSubmatrix<MT,AF,true>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( row_ + m_ == matrix_.rows() ) {
      matrix_.append( row_ + i, column_ + j, value, check );
   }
   else if( !check || !isDefault( value ) ) {
      matrix_.insert( row_ + i, column_ + j, value );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a column.
//
// \param j The index of the column to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a submatrix with elements.
// After completion of column \a j via the append() function, this function can be called to
// finalize column \a j and prepare the next column for insertion process via append().
//
// \note: Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline void SparseSubmatrix<MT,AF,true>::finalize( size_t j )
{
   matrix_.trim( column_ + j );
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
/*!\brief Returns whether the submatrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,AF,true>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,AF,true>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can be used in SMP assignments.
//
// \return \a true in case the submatrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the submatrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline bool SparseSubmatrix<MT,AF,true>::canSMPAssign() const
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,AF,true>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   reserve( 0UL, rows() * columns() );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( IsSymmetric<MT>::value || IsHermitian<MT>::value )
            set( i, j, (~rhs)(i,j) );
         else
            append( i, j, (~rhs)(i,j), true );
      }
      finalize( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the sparse matrix
        , bool AF >       // Alignment flag
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,AF,true>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   reserve( 0UL, (~rhs).nonZeros() );

   for( size_t j=0UL; j<(~rhs).columns(); ++j ) {
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element ) {
         if( IsSymmetric<MT>::value || IsHermitian<MT>::value )
            set( element->index(), j, element->value() );
         else
            append( element->index(), j, element->value(), true );
      }
      finalize( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the sparse matrix
        , bool AF >       // Alignment flag
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,AF,true>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   typedef typename MT2::ConstIterator  RhsIterator;

   // Counting the number of elements per column
   std::vector<size_t> columnLengths( n_, 0UL );
   for( size_t i=0UL; i<m_; ++i ) {
      for( RhsIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         ++columnLengths[element->index()];
   }

   // Resizing the sparse matrix
   for( size_t j=0UL; j<n_; ++j ) {
      reserve( j, columnLengths[j] );
   }

   // Appending the elements to the columns of the sparse matrix
   for( size_t i=0UL; i<m_; ++i ) {
      for( RhsIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         if( IsSymmetric<MT>::value || IsHermitian<MT>::value )
            set( i, element->index(), element->value() );
         else
            append( i, element->index(), element->value(), true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,AF,true>::addAssign( const DenseMatrix<MT2,SO>& rhs )
{
   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const AddType tmp( serial( *this + (~rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO >     // Storage order of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,AF,true>::addAssign( const SparseMatrix<MT2,SO>& rhs )
{
   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const AddType tmp( serial( *this + (~rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,AF,true>::subAssign( const DenseMatrix<MT2,SO>& rhs )
{
   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const SubType tmp( serial( *this - (~rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO >     // Storage order of the right-hand sparse matrix
inline void SparseSubmatrix<MT,AF,true>::subAssign( const SparseMatrix<MT2,SO>& rhs )
{
   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   const SubType tmp( serial( *this - (~rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  SPARSESUBMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseSubmatrix operators */
//@{
template< typename MT, bool AF, bool SO >
inline void reset( SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline void reset( SparseSubmatrix<MT,AF,SO>& sm, size_t i );

template< typename MT, bool AF, bool SO >
inline void clear( SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isDefault( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isIntact( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isSymmetric( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isHermitian( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isLower( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isUniLower( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isStrictlyLower( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isUpper( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isUniUpper( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isStrictlyUpper( const SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isSame( const SparseSubmatrix<MT,AF,SO>& a, const SparseMatrix<MT,SO>& b );

template< typename MT, bool AF, bool SO >
inline bool isSame( const SparseMatrix<MT,SO>& a, const SparseSubmatrix<MT,AF,SO>& b );

template< typename MT, bool AF, bool SO >
inline bool isSame( const SparseSubmatrix<MT,AF,SO>& a, const SparseSubmatrix<MT,AF,SO>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given sparse submatrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be resetted.
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void reset( SparseSubmatrix<MT,AF,SO>& sm )
{
   sm.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given sparse submatrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given sparse submatrix to
// their default value. In case the given submatrix is a \a rowMajor matrix the function resets
// the values in row \a i, if it is a \a columnMajor matrix the function resets the values in
// column \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void reset( SparseSubmatrix<MT,AF,SO>& sm, size_t i )
{
   sm.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given sparse matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse matrix to be cleared.
// \return void
//
// Clearing a sparse submatrix is equivalent to resetting it via the reset() function.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void clear( SparseSubmatrix<MT,AF,SO>& sm )
{
   sm.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse submatrix is in default state.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be tested for its default state.
// \return \a true in case the given submatrix is component-wise zero, \a false otherwise.
//
// This function checks whether the submatrix is in default state. For instance, in
// case the submatrix is instantiated for a built-in integral or floating point data type, the
// function returns \a true in case all submatrix elements are 0 and \a false in case any submatrix
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isDefault( const SparseSubmatrix<MT,AF,SO>& sm )
{
   using blaze::isDefault;

   typedef typename SparseSubmatrix<MT,AF,SO>::ConstIterator  ConstIterator;

   const size_t iend( ( SO == rowMajor)?( sm.rows() ):( sm.columns() ) );

   for( size_t i=0UL; i<iend; ++i ) {
      for( ConstIterator element=sm.begin(i); element!=sm.end(i); ++element )
         if( !isDefault( element->value() ) ) return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given sparse submatrix are intact.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be tested.
// \return \a true in case the given submatrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the sparse submatrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isIntact( const SparseSubmatrix<MT,AF,SO>& sm )
{
   return ( sm.row_ + sm.m_ <= sm.matrix_.rows() &&
            sm.column_ + sm.n_ <= sm.matrix_.columns() &&
            isIntact( sm.matrix_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is symmetric.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is symmetric, \a false if not.
//
// This function checks if the given sparse submatrix is symmetric. The submatrix is considered to
// be symmetric if it is a square matrix whose transpose is equal to itself (\f$ A = A^T \f$). The
// following code example demonstrates the use of the function:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isSymmetric( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSymmetric( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsSymmetric<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isSymmetric( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is Hermitian.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is Hermitian, \a false if not.
//
// This function checks if the given sparse submatrix is Hermitian. The submatrix is considered
// to be Hermitian if it is a square matrix whose transpose is equal to its conjugate transpose
// (\f$ A = \overline{A^T} \f$). The following code example demonstrates the use of the function:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isHermitian( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isHermitian( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsHermitian<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isHermitian( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is a lower triangular matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is a lower triangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is a lower triangular matrix. The matrix is
// considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        l_{1,0} & l_{1,1} & 0       & \cdots & 0       \\
                        l_{2,0} & l_{2,1} & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially lower triangular.
// The following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isLower( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isLower( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsLower<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isLower( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is a lower unitriangular matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is a lower unitriangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is a lower unitriangular matrix. The matrix
// is considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 1       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 1       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 1      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isUniLower( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isUniLower( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsUniLower<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isUniLower( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is a strictly lower triangular matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is a strictly lower triangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is a strictly lower triangular matrix. The
// matrix is considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 0       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 0       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 0      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isStrictlyLower( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isStrictlyLower( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsStrictlyLower<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isStrictlyLower( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is an upper triangular matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is an upper triangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is an upper triangular matrix. The matrix
// is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        u_{0,0} & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & u_{1,1} & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & u_{2,2} & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & u_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially upper triangular.
// The following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isUpper( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsUpper<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isUpper( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is an upper unitriangular matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is an upper unitriangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is an upper triangular matrix. The matrix
// is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 1       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 1       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 1       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isUniUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isUniUpper( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsUniUpper<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isUniUpper( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse submatrix is a strictly upper triangular matrix.
// \ingroup sparse_submatrix
//
// \param sm The sparse submatrix to be checked.
// \return \a true if the submatrix is a strictly upper triangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is a strictly upper triangular matrix. The
// matrix is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 0       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 0       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 0       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::SparseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isStrictlyUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isStrictlyUpper( const SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef SparseMatrix< SparseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsStrictlyUpper<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isStrictlyUpper( static_cast<const BaseType&>( sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse matrix and submatrix represent the same observable state.
// \ingroup sparse_submatrix
//
// \param a The sparse submatrix to be tested for its state.
// \param b The sparse matrix to be tested for its state.
// \return \a true in case the sparse submatrix and matrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given submatrix refers to the full given
// sparse matrix and by that represents the same observable state. In this case, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSame( const SparseSubmatrix<MT,AF,SO>& a, const SparseMatrix<MT,SO>& b )
{
   return ( isSame( a.matrix_, ~b ) && ( a.rows() == (~b).rows() ) && ( a.columns() == (~b).columns() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse matrix and submatrix represent the same observable state.
// \ingroup sparse_submatrix
//
// \param a The sparse matrix to be tested for its state.
// \param b The sparse submatrix to be tested for its state.
// \return \a true in case the sparse matrix and submatrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given submatrix refers to the full given
// sparse matrix and by that represents the same observable state. In this case, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSame( const SparseMatrix<MT,SO>& a, const SparseSubmatrix<MT,AF,SO>& b )
{
   return ( isSame( ~a, b.matrix_ ) && ( (~a).rows() == b.rows() ) && ( (~a).columns() == b.columns() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given submatrices represent the same observable state.
// \ingroup sparse_submatrix
//
// \param a The first sparse submatrix to be tested for its state.
// \param b The second sparse submatrix to be tested for its state.
// \return \a true in case the two submatrices share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given submatrices refer to exactly the
// same part of the same sparse matrix. In case both submatrices represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSame( const SparseSubmatrix<MT,AF,SO>& a, const SparseSubmatrix<MT,AF,SO>& b )
{
   return ( isSame( a.matrix_, b.matrix_ ) &&
            ( a.row_ == b.row_ ) && ( a.column_ == b.column_ ) &&
            ( a.m_ == b.m_ ) && ( a.n_ == b.n_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a sparse submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO      // Storage order
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryAssign( const SparseSubmatrix<MT,AF,SO>& lhs, const Vector<VT,TF>& rhs,
                       size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a matrix to a sparse submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO1      // Storage order
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryAssign( const SparseSubmatrix<MT1,AF,SO1>& lhs, const Matrix<MT2,SO2>& rhs,
                       size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return tryAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a sparse submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO      // Storage order
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const SparseSubmatrix<MT,AF,SO>& lhs, const Vector<VT,TF>& rhs,
                          size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryAddAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a sparse submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side matrix to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO1      // Storage order
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryAddAssign( const SparseSubmatrix<MT1,AF,SO1>& lhs, const Matrix<MT2,SO2>& rhs,
                          size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return tryAddAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a sparse
//        submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO      // Storage order
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool trySubAssign( const SparseSubmatrix<MT,AF,SO>& lhs, const Vector<VT,TF>& rhs,
                          size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return trySubAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a sparse
//        submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side matrix to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO1      // Storage order
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool trySubAssign( const SparseSubmatrix<MT1,AF,SO1>& lhs, const Matrix<MT2,SO2>& rhs,
                          size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return trySubAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a sparse
//        submatrix.
// \ingroup sparse_submatrix
//
// \param lhs The target left-hand side sparse submatrix.
// \param rhs The right-hand side vector to be multiplied.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO      // Storage order
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const SparseSubmatrix<MT,AF,SO>& lhs, const Vector<VT,TF>& rhs,
                           size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryMultAssign( lhs.matrix_, ~rhs, lhs.row_ + row, lhs.column_ + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given sparse submatrix.
// \ingroup sparse_submatrix
//
// \param sm The submatrix to be derestricted.
// \return Submatrix without access restrictions.
//
// This function removes all restrictions on the data access to the given submatrix. It returns a
// submatrix that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DerestrictTrait< SparseSubmatrix<MT,AF,SO> >::Type
   derestrict( SparseSubmatrix<MT,AF,SO>& sm )
{
   typedef typename DerestrictTrait< SparseSubmatrix<MT,AF,SO> >::Type  ReturnType;
   return ReturnType( derestrict( sm.matrix_ ), sm.row_, sm.column_, sm.m_, sm.n_ );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another sparse submatrix.
// \ingroup views
//
// \param sm The constant sparse submatrix
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the other sparse submatrix.
//
// This function returns an expression representing the specified submatrix of the given
// sparse submatrix.
*/
template< bool AF1     // Required alignment flag
        , typename MT  // Type of the sparse submatrix
        , bool AF2     // Present alignment flag
        , bool SO >    // Storage order
inline const SparseSubmatrix<MT,AF1,SO>
   submatrix( const SparseSubmatrix<MT,AF2,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   if( ( row + m > sm.rows() ) || ( column + n > sm.columns() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
   }

   return SparseSubmatrix<MT,AF1,SO>( sm.matrix_, sm.row_ + row, sm.column_ + column, m, n );
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
template< typename MT, bool AF, bool SO >
struct IsRestricted< SparseSubmatrix<MT,AF,SO> > : public If< IsRestricted<MT>, TrueType, FalseType >::Type
{
   enum { value = IsRestricted<MT>::value };
   typedef typename If< IsRestricted<MT>, TrueType, FalseType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DERESTRICTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF, bool SO >
struct DerestrictTrait< SparseSubmatrix<MT,AF,SO> >
{
   typedef SparseSubmatrix< typename RemoveReference< typename DerestrictTrait<MT>::Type >::Type, AF >  Type;
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
template< typename MT, bool AF, bool SO, typename T >
struct AddTrait< SparseSubmatrix<MT,AF,SO>, T >
{
   typedef typename AddTrait< typename SubmatrixTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool AF, bool TF >
struct AddTrait< T, SparseSubmatrix<MT,AF,TF> >
{
   typedef typename AddTrait< T, typename SubmatrixTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool AF, bool SO, typename T >
struct SubTrait< SparseSubmatrix<MT,AF,SO>, T >
{
   typedef typename SubTrait< typename SubmatrixTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool AF, bool TF >
struct SubTrait< T, SparseSubmatrix<MT,AF,TF> >
{
   typedef typename SubTrait< T, typename SubmatrixTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool AF, bool SO, typename T >
struct MultTrait< SparseSubmatrix<MT,AF,SO>, T >
{
   typedef typename MultTrait< typename SubmatrixTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool AF, bool TF >
struct MultTrait< T, SparseSubmatrix<MT,AF,TF> >
{
   typedef typename MultTrait< T, typename SubmatrixTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool AF, bool SO, typename T >
struct DivTrait< SparseSubmatrix<MT,AF,SO>, T >
{
   typedef typename DivTrait< typename SubmatrixTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool AF, bool TF >
struct DivTrait< T, SparseSubmatrix<MT,AF,TF> >
{
   typedef typename DivTrait< T, typename SubmatrixTrait<MT>::Type >::Type  Type;
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
template< typename MT, bool AF, bool SO >
struct SubmatrixTrait< SparseSubmatrix<MT,AF,SO> >
{
   typedef typename SubmatrixTrait< typename SparseSubmatrix<MT,AF,SO>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXEXPRTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< SparseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef SparseSubmatrix<MT,AF2,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< const SparseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef SparseSubmatrix<MT,AF2,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< volatile SparseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef SparseSubmatrix<MT,AF2,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< const volatile SparseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef SparseSubmatrix<MT,AF2,SO>  Type;
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
template< typename MT, bool AF, bool SO >
struct RowTrait< SparseSubmatrix<MT,AF,SO> >
{
   typedef typename RowTrait< typename SparseSubmatrix<MT,AF,SO>::ResultType >::Type  Type;
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
template< typename MT, bool AF, bool SO >
struct ColumnTrait< SparseSubmatrix<MT,AF,SO> >
{
   typedef typename ColumnTrait< typename SparseSubmatrix<MT,AF,SO>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
