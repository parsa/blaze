//=================================================================================================
/*!
//  \file blaze/math/views/DenseSubmatrix.h
//  \brief Header file for the DenseSubmatrix class template
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

#ifndef _BLAZE_MATH_VIEWS_DENSESUBMATRIX_H_
#define _BLAZE_MATH_VIEWS_DENSESUBMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <stdexcept>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Restricted.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Submatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/AlignmentFlag.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Streaming.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/AlignedArray.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_submatrix Dense Submatrix
// \ingroup views
*/
/*!\brief View on a specific submatrix of a dense matrix.
// \ingroup dense_submatrix
//
// The DenseSubmatrix template represents a view on a specific submatrix of a dense matrix
// primitive. The type of the dense matrix is specified via the first template parameter:

   \code
   template< typename MT, bool AF, bool SO >
   class DenseSubmatrix;
   \endcode

//  - MT: specifies the type of the dense matrix primitive. DenseSubmatrix can be used with every
//        dense matrix primitive, but does not work with any matrix expression type.
//  - AF: the alignment flag specifies whether the submatrix is aligned (\a blaze::aligned) or
//        unaligned (\a blaze::unaligned). The default value is \a blaze::unaligned.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the dense matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//
//
// \n \section dense_submatrix_setup Setup of Dense Submatrices
//
// A view on a dense submatrix can be created very conveniently via the \c submatrix() function:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrixType;

   DenseMatrixType A;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x16, starting in row 0 and column 4
   blaze::DenseSubmatrix<DenseMatrixType> sm = submatrix( A, 0UL, 4UL, 8UL, 16UL );
   \endcode

// This view can be treated as any other dense matrix, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The view can also be used on both
// sides of an assignment: The submatrix can either be used as an alias to grant write access to
// a specific submatrix of a dense matrix primitive on the left-hand side of an assignment or
// to grant read-access to a specific submatrix of a dense matrix primitive or expression on
// the right-hand side of an assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicMatrix<double,blaze::columnMajor>  DenseMatrixType;
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;

   DenseMatrixType A, B;
   SparseMatrixType C;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 2
   blaze::DenseSubmatrix<DenseMatrixType> sm = submatrix( A, 0UL, 2UL, 8UL, 4UL );

   // Setting the submatrix of A to a 8x4 submatrix of B
   sm = submatrix( B, 0UL, 0UL, 8UL, 4UL );

   // Copying the sparse matrix C into another 8x4 submatrix of A
   submatrix( A, 8UL, 2UL, 8UL, 4UL ) = C;

   // Assigning part of the result of a matrix addition to the first submatrix
   sm = submatrix( B + C, 0UL, 0UL, 8UL, 4UL );
   \endcode

// \n \section dense_submatrix_element_access Element access
//
// A dense submatrix can be used like any other dense matrix. For instance, the elements of the
// dense submatrix can be directly accessed with the function call operator:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a 8x8 submatrix, starting from position (4,4)
   blaze::DenseSubmatrix<MatrixType> sm = submatrix( A, 4UL, 4UL, 8UL, 8UL );

   // Setting the element (0,0) of the submatrix, which corresponds to
   // the element at position (4,4) in matrix A
   sm(0,0) = 2.0;
   \endcode

// Alternatively, the elements of a submatrix can be traversed via (const) iterators. Just as
// with matrices, in case of non-const submatrices, \c begin() and \c end() return an Iterator,
// which allows a manipulation of the non-zero values, in case of constant submatrices a
// ConstIterator is returned:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>          SubmatrixType;

   MatrixType A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a specific submatrix of matrix A
   SubmatrixType sm = submatrix( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( SubmatrixType::Iterator it=sm.begin(0); it!=sm.end(0); ++it ) {
      *it = ...;  // OK: Write access to the dense submatrix value.
      ... = *it;  // OK: Read access to the dense submatrix value.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( SubmatrixType::ConstIterator it=sm.begin(1); it!=sm.end(1); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense submatrix value.
   }
   \endcode

// \n \section dense_submatrix_common_operations Common Operations
//
// The current size of the matrix, i.e. the number of rows or columns can be obtained via the
// \c rows() and \c columns() functions, the current total capacity via the \c capacity() function,
// and the number of non-zero elements via the \c nonZeros() function. However, since submatrices
// are views on a specific submatrix of a matrix, several operations are not possible on views,
// such as resizing and swapping:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>          SubmatrixType;

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

// \n \section dense_submatrix_arithmetic_operations Arithmetic Operations
//
// The following example gives an impression of the use of DenseSubmatrix within arithmetic
// operations. All operations (addition, subtraction, multiplication, scaling, ...) can be
// performed on all possible combinations of dense and sparse matrices with fitting element
// types:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>     DenseMatrixType;
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;
   DenseMatrixType D1, D2, D3;
   SparseMatrixType S1, S2;

   typedef blaze::CompressedVector<double,blaze::columnVector>  SparseVectorType;
   SparseVectorType a, b;

   // ... Resizing and initialization

   typedef DenseSubmatrix<DenseMatrixType>  SubmatrixType;
   SubmatrixType sm = submatrix( D1, 0UL, 0UL, 8UL, 8UL );  // View on the 8x8 submatrix of matrix D1
                                                            // starting from row 0 and column 0

   submatrix( D1, 0UL, 8UL, 8UL, 8UL ) = D2;  // Dense matrix initialization of the 8x8 submatrix
                                              // starting in row 0 and column 8
   sm = S1;                                   // Sparse matrix initialization of the second 8x8 submatrix

   D3 = sm + D2;                                    // Dense matrix/dense matrix addition
   S2 = S1  - submatrix( D1, 8UL, 0UL, 8UL, 8UL );  // Sparse matrix/dense matrix subtraction
   D2 = sm * submatrix( D1, 8UL, 8UL, 8UL, 8UL );   // Dense matrix/dense matrix multiplication

   submatrix( D1, 8UL, 0UL, 8UL, 8UL ) *= 2.0;      // In-place scaling of a submatrix of D1
   D2 = submatrix( D1, 8UL, 8UL, 8UL, 8UL ) * 2.0;  // Scaling of the a submatrix of D1
   D2 = 2.0 * sm;                                   // Scaling of the a submatrix of D1

   submatrix( D1, 0UL, 8UL, 8UL, 8UL ) += D2;  // Addition assignment
   submatrix( D1, 8UL, 0UL, 8UL, 8UL ) -= S1;  // Subtraction assignment
   submatrix( D1, 8UL, 8UL, 8UL, 8UL ) *= sm;  // Multiplication assignment

   a = submatrix( D1, 4UL, 4UL, 8UL, 8UL ) * b;  // Dense matrix/sparse vector multiplication
   \endcode

// \n \section dense_submatrix_aligned_submatrix Aligned Submatrices
//
// Usually submatrices can be defined anywhere within a matrix. They may start at any position and
// may have an arbitrary extension (only restricted by the extension of the underlying matrix).
// However, in contrast to matrices themselves, which are always properly aligned in memory and
// therefore can provide maximum performance, this means that submatrices in general have to be
// considered to be unaligned. This can be made explicit by the \a blaze::unaligned flag:

   \code
   using blaze::unaligned;

   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrixType;

   DenseMatrixType A;
   // ... Resizing and initialization

   // Identical creations of an unaligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::DenseSubmatrix<DenseMatrixType>           sm1 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   blaze::DenseSubmatrix<DenseMatrixType>           sm2 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   blaze::DenseSubmatrix<DenseMatrixType,unaligned> sm3 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   blaze::DenseSubmatrix<DenseMatrixType,unaligned> sm4 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// All of these calls to the \c submatrix() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned submatrix. Whereas this may provide
// full flexibility in the creation of submatrices, this might result in performance restrictions
// (even in case the specified submatrix could be aligned). However, it is also possible to create
// aligned submatrices. Aligned submatrices are identical to unaligned submatrices in all aspects,
// except that they may pose additional alignment restrictions and therefore have less flexibility
// during creation, but don't suffer from performance penalties and provide the same performance
// as the underlying matrix. Aligned submatrices are created by explicitly specifying the
// \a blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::DenseSubmatrix<DenseMatrixType,aligned> sv = submatrix<aligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). The following source code gives some
// examples for a double precision dense matrix, assuming that AVX is available, which packs 4
// \c double values into an intrinsic vector:

   \code
   using blaze::rowMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>      MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType,aligned>  SubmatrixType;

   MatrixType D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0) and the number of rows and columns are a multiple of 4
   SubmatrixType dsm1 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 12UL );

   // OK: First row and column and the number of rows and columns are all a multiple of 4
   SubmatrixType dsm2 = submatrix<aligned>( D, 4UL, 12UL, 8UL, 16UL );

   // OK: First row and column are a multiple of 4 and the submatrix includes the last row and column
   SubmatrixType dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First row is not a multiple of 4
   SubmatrixType dsm4 = submatrix<aligned>( D, 2UL, 4UL, 12UL, 12UL );

   // Error: First column is not a multiple of 4
   SubmatrixType dsm5 = submatrix<aligned>( D, 0UL, 2UL, 8UL, 8UL );

   // Error: The number of rows is not a multiple of 4 and the submatrix does not include the last row
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 8UL );

   // Error: The number of columns is not a multiple of 4 and the submatrix does not include the last column
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 11UL );
   \endcode

// \n \section dense_submatrix_on_dense_submatrix Submatrix on Submatrix
//
// It is also possible to create a submatrix view on another submatrix. In this context it is
// important to remember that the type returned by the \c submatrix() function is the same type
// as the type of the given submatrix, since the view on a submatrix is just another view on the
// underlying dense matrix:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>             SubmatrixType;

   MatrixType D1;

   // ... Resizing and initialization

   // Creating a submatrix view on the dense matrix D1
   SubmatrixType sm1 = submatrix( D1, 4UL, 4UL, 8UL, 16UL );

   // Creating a submatrix view on the dense submatrix sm1
   SubmatrixType sm2 = submatrix( sm1, 1UL, 1UL, 4UL, 8UL );
   \endcode

// \n \section dense_submatrix_on_symmetric_matrices Submatrix on Symmetric Matrices
//
// Submatrices can also be created on symmetric matrices (see the SymmetricMatrix class template):

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::DenseSubmatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >   SymmetricDynamicType;
   typedef DenseSubmatrix< SymmetricDynamicType >  SubmatrixType;

   // Setup of a 16x16 symmetric matrix
   SymmetricDynamicType A( 16UL );

   // Creating a dense submatrix of size 8x12, starting in row 2 and column 4
   SubmatrixType sm = submatrix( A, 2UL, 4UL, 8UL, 12UL );
   \endcode

// It is important to note, however, that (compound) assignments to such submatrices have a
// special restriction: The symmetry of the underlying symmetric matrix must not be broken!
// Since the modification of element \f$ a_{ij} \f$ of a symmetric matrix also modifies the
// element \f$ a_{ji} \f$, the matrix to be assigned must be structured such that the symmetry
// of the symmetric matrix is preserved. Otherwise a \a std::invalid_argument exception is
// thrown:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Setup of two default 4x4 symmetric matrices
   SymmetricMatrix< DynamicMatrix<int> > A1( 4 ), A2( 4 );

   // Setup of the 3x2 dynamic matrix
   //
   //       ( 0 9 )
   //   B = ( 9 8 )
   //       ( 0 7 )
   //
   DynamicMatrix<int> B( 3UL, 2UL );
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
template< typename MT                                 // Type of the dense matrix
        , bool AF = unaligned                         // Alignment flag
        , bool SO = IsColumnMajorMatrix<MT>::value >  // Storage order
class DenseSubmatrix : public DenseMatrix< DenseSubmatrix<MT,AF,SO>, SO >
                     , private Submatrix
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the matrix element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense submatrix will
       return references and iterators to const. Otherwise \a useConst will be set to 0 and
       the dense submatrix will offer write access to the dense matrix elements both via
       the function call operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubmatrix<MT,AF,SO>            This;           //!< Type of this DenseSubmatrix instance.
   typedef typename SubmatrixTrait<MT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the submatrix elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the submatrix elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubmatrix&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant submatrix value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstPointer, ElementType* >::Type  Pointer;
   //**********************************************************************************************

   //**SubmatrixIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse submatrix.
   */
   template< typename IteratorType >  // Type of the dense matrix iterator
   class SubmatrixIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! The iterator category.
      typedef typename std::iterator_traits<IteratorType>::iterator_category  IteratorCategory;

      //! Type of the underlying elements.
      typedef typename std::iterator_traits<IteratorType>::value_type  ValueType;

      //! Pointer return type.
      typedef typename std::iterator_traits<IteratorType>::pointer  PointerType;

      //! Reference return type.
      typedef typename std::iterator_traits<IteratorType>::reference  ReferenceType;

      //! Difference between two iterators.
      typedef typename std::iterator_traits<IteratorType>::difference_type  DifferenceType;

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the SubmatrixIterator class.
      */
      inline SubmatrixIterator()
         : iterator_ (       )  // Iterator to the current submatrix element
         , final_    ( 0UL   )  // The final iterator for intrinsic operations
         , rest_     ( 0UL   )  // The number of remaining elements beyond the final iterator
         , isAligned_( false )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the SubmatrixIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param finalIterator The final iterator for intrinsic operations.
      // \param remainingElements The number of remaining elements beyond the final iterator.
      // \param isMemoryAligned Memory alignment flag.
      */
      inline SubmatrixIterator( IteratorType iterator, IteratorType finalIterator,
                                size_t remainingElements, bool isMemoryAligned )
         : iterator_ ( iterator          )  // Iterator to the current submatrix element
         , final_    ( finalIterator     )  // The final iterator for intrinsic operations
         , rest_     ( remainingElements )  // The number of remaining elements beyond the final iterator
         , isAligned_( isMemoryAligned   )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename IteratorType2 >
      inline SubmatrixIterator( const SubmatrixIterator<IteratorType2>& it )
         : iterator_ ( it.base()      )  // Iterator to the current submatrix element
         , final_    ( it.final()     )  // The final iterator for intrinsic operations
         , rest_     ( it.rest()      )  // The number of remaining elements beyond the final iterator
         , isAligned_( it.isAligned() )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline SubmatrixIterator& operator+=( size_t inc ) {
         iterator_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline SubmatrixIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubmatrixIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator++( int ) {
         return SubmatrixIterator( iterator_++, final_, rest_, isAligned_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline SubmatrixIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator--( int ) {
         return SubmatrixIterator( iterator_--, final_, rest_, isAligned_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReferenceType operator*() const {
         return *iterator_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Aligned load of an intrinsic element of the dense submatrix.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an aligned load of the current intrinsic element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType load() const {
         return loadu();
      }
      //*******************************************************************************************

      //**Loadu function***************************************************************************
      /*!\brief Unaligned load of an intrinsic element of the dense submatrix.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an unaligned load of the current intrinsic element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType loadu() const {
         if( isAligned_ ) {
            return iterator_.load();
         }
         else if( iterator_ != final_ ) {
            return iterator_.loadu();
         }
         else {
            AlignedArray<ElementType,IT::size> array;
            for( size_t j=0UL; j<rest_; ++j )
               array[j] = *(iterator_+j);
            for( size_t j=rest_; j<IT::size; ++j )
               array[j] = ElementType();
            return blaze::load( array.data() );
         }
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const SubmatrixIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const SubmatrixIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const SubmatrixIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const SubmatrixIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const SubmatrixIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const SubmatrixIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const SubmatrixIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( const SubmatrixIterator& it, size_t inc ) {
         return SubmatrixIterator( it.iterator_ + inc, it.final_, it.rest_, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a SubmatrixIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( size_t inc, const SubmatrixIterator& it ) {
         return SubmatrixIterator( it.iterator_ + inc, it.final_, it.rest_, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const SubmatrixIterator operator-( const SubmatrixIterator& it, size_t dec ) {
         return SubmatrixIterator( it.iterator_ - dec, it.final_, it.rest_, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the submatrix iterator.
      //
      // \return The current position of the submatrix iterator.
      */
      inline IteratorType base() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**Final function***************************************************************************
      /*!\brief Access to the final position of the submatrix iterator.
      //
      // \return The final position of the submatrix iterator.
      */
      inline IteratorType final() const {
         return final_;
      }
      //*******************************************************************************************

      //**Rest function****************************************************************************
      /*!\brief Access to the number of remaining elements beyond the final iterator.
      //
      // \return The number of remaining elements beyond the final iterator.
      */
      inline size_t rest() const {
         return rest_;
      }
      //*******************************************************************************************

      //**IsAligned function***********************************************************************
      /*!\brief Access to the iterator's memory alignment flag.
      //
      // \return \a true in case the iterator is aligned, \a false if it is not.
      */
      inline bool isAligned() const {
         return isAligned_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;   //!< Iterator to the current submatrix element.
      IteratorType final_;      //!< The final iterator for intrinsic operations.
      size_t       rest_;       //!< The number of remaining elements beyond the final iterator.
      bool         isAligned_;  //!< Memory alignment flag.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubmatrixIterator<typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename IfTrue< useConst, ConstIterator, SubmatrixIterator<typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = MT::vectorizable };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseSubmatrix& operator=( const ElementType& rhs );
   inline DenseSubmatrix& operator=( const DenseSubmatrix& rhs );

   template< typename MT2, bool SO2 >
   inline DenseSubmatrix& operator=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline DenseSubmatrix& operator*=( const Matrix<MT2,SO2>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          rows() const;
                              inline size_t          columns() const;
                              inline size_t          spacing() const;
                              inline size_t          capacity() const;
                              inline size_t          capacity( size_t i ) const;
                              inline size_t          nonZeros() const;
                              inline size_t          nonZeros( size_t i ) const;
                              inline void            reset();
                              inline void            reset( size_t i );
                              inline DenseSubmatrix& transpose();
   template< typename Other > inline DenseSubmatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load ( size_t i, size_t j ) const;
   inline IntrinsicType loadu( size_t i, size_t j ) const;

   inline void store ( size_t i, size_t j, const IntrinsicType& value );
   inline void storeu( size_t i, size_t j, const IntrinsicType& value );
   inline void stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,SO>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,SO>& rhs );

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,!SO>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,SO>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,!SO>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,SO>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,SO>& rhs );

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,!SO>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,SO>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,!SO>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,SO>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,SO>& rhs );

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,!SO>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,SO>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,!SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const;

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;   //!< The dense matrix containing the submatrix.
   const size_t row_;      //!< The first row of the submatrix.
   const size_t column_;   //!< The first column of the submatrix.
   const size_t m_;        //!< The number of rows of the submatrix.
   const size_t n_;        //!< The number of columns of the submatrix.
   const size_t rest_;     //!< The number of remaining elements in an unaligned intrinsic operation.
   const size_t final_;    //!< The final index for unaligned intrinsic operations.
                           /*!< In case the submatrix is not fully aligned and the submatrix is
                                involved in a vectorized operation, the final index indicates at
                                which index a special treatment for the remaining elements is
                                required. */
   const bool isAligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the submatrix is fully aligned.
                                In case the submatrix is fully aligned, no special handling has to
                                be used for the last elements of the submatrix in a vectorized
                                operation. In order to be aligned, the following conditions must
                                hold for the submatrix:
                                 - The first element of each row/column must be aligned
                                 - The submatrix must be at the end of the given matrix or
                                 - The number of rows/columns of the submatrix must be a multiple
                                   of the number of values per intrinsic element. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename MT2, bool AF2, bool SO2 > friend class DenseSubmatrix;

   template< bool AF1, typename MT2, bool AF2, bool SO2 >
   friend const DenseSubmatrix<MT2,AF1,SO2>
      submatrix( const DenseSubmatrix<MT2,AF2,SO2>& dm, size_t row, size_t column, size_t m, size_t n );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSymmetric( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isLower( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isUpper( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseMatrix<MT2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseMatrix<MT2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend typename DerestrictTrait< DenseSubmatrix<MT2,AF2,SO2> >::Type
      derestrict( DenseSubmatrix<MT2,AF2,SO2>& dm );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
   BLAZE_STATIC_ASSERT( !IsRestricted<MT>::value || IsLower<MT>::value || IsUpper<MT>::value );
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
/*!\brief The constructor for DenseSubmatrix.
//
// \param matrix The dense matrix containing the submatrix.
// \param row The index of the first row of the submatrix in the given dense matrix.
// \param column The index of the first column of the submatrix in the given dense matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given dense matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline DenseSubmatrix<MT,AF,SO>::DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_   ( matrix       )  // The dense matrix containing the submatrix
   , row_      ( row          )  // The first row of the submatrix
   , column_   ( column       )  // The first column of the submatrix
   , m_        ( m            )  // The number of rows of the submatrix
   , n_        ( n            )  // The number of columns of the submatrix
   , rest_     ( n % IT::size )  // The number of remaining elements in an unaligned intrinsic operation
   , final_    ( n - rest_    )  // The final index for unaligned intrinsic operations
   , isAligned_( ( column % IT::size == 0UL ) &&
                 ( column + n == matrix.columns() || n % IT::size == 0UL ) )
{
   if( ( row + m > matrix.rows() ) || ( column + n > matrix.columns() ) )
      throw std::invalid_argument( "Invalid submatrix specification" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::Reference
   DenseSubmatrix<MT,AF,SO>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row_+i,column_+j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::ConstReference
   DenseSubmatrix<MT,AF,SO>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row_+i,column_+j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::Pointer DenseSubmatrix<MT,AF,SO>::data()
{
   return matrix_.data() + row_*spacing() + column_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::ConstPointer DenseSubmatrix<MT,AF,SO>::data() const
{
   return matrix_.data() + row_*spacing() + column_;
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::Iterator DenseSubmatrix<MT,AF,SO>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   const typename MT::Iterator first( matrix_.begin( row_ + i ) + column_ );
   return Iterator( first, first + final_, rest_, isAligned_ );
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::ConstIterator
   DenseSubmatrix<MT,AF,SO>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   const typename MT::ConstIterator first( matrix_.cbegin( row_ + i ) + column_ );
   return ConstIterator( first, first + final_, rest_, isAligned_ );
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::ConstIterator
   DenseSubmatrix<MT,AF,SO>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   const typename MT::ConstIterator first( matrix_.cbegin( row_ + i ) + column_ );
   return ConstIterator( first, first + final_, rest_, isAligned_ );
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::Iterator DenseSubmatrix<MT,AF,SO>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   const typename MT::Iterator last( matrix_.begin( row_ + i ) + column_ + n_ );
   return Iterator( last, last, rest_, isAligned_ );
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::ConstIterator
   DenseSubmatrix<MT,AF,SO>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   const typename MT::ConstIterator last( matrix_.cbegin( row_ + i ) + column_ + n_ );
   return ConstIterator( last, last, rest_, isAligned_ );
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::ConstIterator
   DenseSubmatrix<MT,AF,SO>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   const typename MT::ConstIterator last( matrix_.cbegin( row_ + i ) + column_ + n_ );
   return ConstIterator( last, last, rest_, isAligned_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline DenseSubmatrix<MT,AF,SO>& DenseSubmatrix<MT,AF,SO>::operator=( const ElementType& rhs )
{
   const size_t iend( row_ + m_ );

   for( size_t i=row_; i<iend; ++i )
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

      for( size_t j=jbegin; j<jend; ++j )
         matrix_(i,j) = rhs;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DenseSubmatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline DenseSubmatrix<MT,AF,SO>& DenseSubmatrix<MT,AF,SO>::operator=( const DenseSubmatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ && column_ == rhs.column_ ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() )
      throw std::invalid_argument( "Submatrix sizes do not match" );

   if( !preservesInvariant( matrix_, rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,AF,SO>& DenseSubmatrix<MT,AF,SO>::operator=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   typedef typename If< IsAdaptor<MT>, typename MT2::CompositeType, const MT2& >::Type  Right;
   Right right( ~rhs );

   if( !preservesInvariant( matrix_, right ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   if( IsSparseMatrix<MT2>::value )
      reset();

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,AF,SO>& >::Type
   DenseSubmatrix<MT,AF,SO>::operator+=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const AddType tmp( *this + (~rhs) );
      smpAssign( left, tmp );
   }
   else {
      smpAddAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,AF,SO>& >::Type
   DenseSubmatrix<MT,AF,SO>::operator+=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const AddType tmp( *this + (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,AF,SO>& >::Type
   DenseSubmatrix<MT,AF,SO>::operator-=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const SubType tmp( *this - (~rhs ) );
      smpAssign( left, tmp );
   }
   else {
      smpSubAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,AF,SO>& >::Type
   DenseSubmatrix<MT,AF,SO>::operator-=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const SubType tmp( *this - (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,AF,SO>& DenseSubmatrix<MT,AF,SO>::operator*=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const MultType tmp( *this * (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a dense submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
*/
template< typename MT       // Type of the dense matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,AF,SO> >::Type&
   DenseSubmatrix<MT,AF,SO>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) * rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense submatrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the dense matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,AF,SO> >::Type&
   DenseSubmatrix<MT,AF,SO>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) / rhs );

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the number of rows of the dense submatrix.
//
// \return The number of rows of the dense submatrix.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::rows() const
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of columns of the dense submatrix.
//
// \return The number of columns of the dense submatrix.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::columns() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows/columns.
//
// \return The spacing between the beginning of two rows/columns.
//
// This function returns the spacing between the beginning of two rows/columns, i.e. the
// total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::spacing() const
{
   return matrix_.spacing();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::capacity() const
{
   return rows() * columns();
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::capacity( size_t i ) const
{
   UNUSED_PARAMETER( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::nonZeros() const
{
   const size_t iend( row_ + m_ );
   const size_t jend( column_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t i=row_; i<iend; ++i )
      for( size_t j=column_; j<jend; ++j )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline size_t DenseSubmatrix<MT,AF,SO>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( column_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t j=column_; j<jend; ++j )
      if( !isDefault( matrix_(row_+i,j) ) )
         ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void DenseSubmatrix<MT,AF,SO>::reset()
{
   using blaze::clear;

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

      for( size_t j=jbegin; j<jend; ++j )
         clear( matrix_(i,j) );
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void DenseSubmatrix<MT,AF,SO>::reset( size_t i )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

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

   for( size_t j=jbegin; j<jend; ++j )
      clear( matrix_(row_+i,j) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transposing the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::runtime_error Invalid transpose of a non-quadratic submatrix.
// \exception std::runtime_error Invalid transpose of a lower matrix.
// \exception std::runtime_error Invalid transpose of an upper matrix.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// in case the underlying matrix is a lower or upper triangular matrix the function can only be
// used in case the submatrix does not contain elements from the upper or lower part of the matrix,
// respectively. The attempt to transpose a non-quadratic submatrix or an invalid part of a lower
// or triangular matrix results in a \a std::runtime_error exception.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline DenseSubmatrix<MT,AF,SO>& DenseSubmatrix<MT,AF,SO>::transpose()
{
   if( rows() != columns() )
      throw std::runtime_error( "Invalid transpose of a non-quadratic submatrix" );

   if( IsLower<MT>::value && ( row_ + 1UL < column_ + n_ ) )
      throw std::runtime_error( "Invalid transpose of a lower matrix" );

   if( IsUpper<MT>::value && ( column_ + 1UL < row_ + m_ ) )
      throw std::runtime_error( "Invalid transpose of an upper matrix" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( trans(*this) );
   smpAssign( left, tmp );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales all elements of the submatrix by the given scalar value \a scalar. Note
// that the function cannot be used to scale a submatrix on a lower or upper unitriangular matrix.
// The attempt to scale such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the dense matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline DenseSubmatrix<MT,AF,SO>& DenseSubmatrix<MT,AF,SO>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t iend( row_ + m_ );

   for( size_t i=row_; i<iend; ++i )
   {
      const size_t jbegin( ( IsUpper<MT>::value )
                           ?( ( IsStrictlyUpper<MT>::value )
                              ?( max( i+1UL, column_ ) )
                              :( max( i, column_ ) ) )
                           :( column_ ) );
      const size_t jend  ( ( IsLower<MT>::value )
                           ?( ( IsStrictlyLower<MT>::value )
                              ?( min( i, column_+n_ ) )
                              :( min( i+1UL, column_+n_ ) ) )
                           :( column_+n_ ) );

      for( size_t j=jbegin; j<jend; ++j )
         matrix_(i,j) *= scalar;
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
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool DenseSubmatrix<MT,AF,SO>::hasOverlap() const
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric<MT>::value, "Unsymmetric matrix detected" );

   if( ( row_ + m_ <= column_ ) || ( column_ + n_ <= row_ ) )
      return false;
   else return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying matrix.
//
// \param lhs The matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying matrix of type \a MT would be
// violated by an assignment of the given matrix \a rhs. In case the matrix would be preserved,
// the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the left-hand side dense matrix
        , bool SO2      // Storage order of the left-hand side dense matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3 >    // Storage order of the right-hand side matrix
inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   UNUSED_PARAMETER( lhs, rhs );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying symmetric matrix.
//
// \param lhs The symmetric matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying symmetric matrix of type \a MT would
// be violated by an assignment of the given row-major matrix \a rhs. In case the matrix would be
// preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the left-hand side dense matrix
        , bool SO2      // Storage order of the left-hand side dense matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3 >    // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( !hasOverlap() )
      return true;

   const bool   lower( row_ > column_ );
   const size_t size ( min( row_ + m_, column_ + n_ ) - ( lower ? row_ : column_ ) );

   if( size < 2UL )
      return true;

   const size_t row   ( lower ? 0UL : column_ - row_ );
   const size_t column( lower ? row_ - column_ : 0UL );

   return isSymmetric( submatrix( ~rhs, row, column, size, size ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;

      const size_t jbegin( ( containsDiagonal )
                           ?( ( IsStrictlyLower<MT2>::value )
                              ?( row_ + i - column_ )
                              :( row_ + i - column_ + 1UL ) )
                           :( 0UL ) );

      for( size_t j=jbegin; j<n_; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t ipos( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                         ?( column_ + j - row_ + 1UL )
                         :( column_ + j - row_ ) );
      const size_t iend( min( ipos, m_ ) );

      for( size_t i=0UL; i<iend; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                             ?( row_ + i - column_ )
                             :( row_ + i - column_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(i) );
      RhsIterator element( (~rhs).lowerBound( i, index ) );

      if( IsUniLower<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != row_+i-column_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t index( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                          ?( column_ + j - row_ + 1UL )
                          :( column_ + j - row_ ) );
      const RhsIterator last( (~rhs).lowerBound( min( index, m_ ), j ) );

      if( IsUniLower<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(j) || ( last->index() != column_+j-row_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(j); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t jpos( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                         ?( row_ + i - column_ + 1UL )
                         :( row_ + i - column_ ) );
      const size_t jend( min( jpos, n_ ) );

      for( size_t j=0UL; j<jend; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;

      const size_t ibegin( ( containsDiagonal )
                           ?( ( IsStrictlyUpper<MT2>::value )
                              ?( column_ + j - row_ )
                              :( column_ + j - row_ + 1UL ) )
                           :( 0UL ) );

      for( size_t i=ibegin; i<m_; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t index( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                          ?( row_ + i - column_ + 1UL )
                          :( row_ + i - column_ ) );
      const RhsIterator last( (~rhs).lowerBound( i, min( index, n_ ) ) );

      if( IsUniUpper<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(i) || ( last->index() != row_+i-column_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(i); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                             ?( column_ + j - row_ )
                             :( column_ + j - row_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(j) );
      RhsIterator element( (~rhs).lowerBound( index, j ) );

      if( IsUniUpper<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != column_+j-row_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<n_; ++j ) {
         if( ( row_ + i != column_ + j ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<m_; ++i ) {
         if( ( column_ + j != row_ + i ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t i=0UL; i<m_; ++i ) {
      for( RhsIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element ) {
         if( ( row_ + i != column_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,AF,SO>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t j=0UL; j<n_; ++j ) {
      for( RhsIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element ) {
         if( ( column_ + j != row_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
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
template< typename MT       // Type of the dense matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,AF,SO>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Data type of the foreign dense submatrix
        , bool AF2      // Alignment flag of the foreign dense submatrix
        , bool SO2 >    // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,AF,SO>::canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
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
template< typename MT       // Type of the dense matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,AF,SO>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT   // Type of the dense matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Data type of the foreign dense submatrix
        , bool AF2      // Alignment flag of the foreign dense submatrix
        , bool SO2 >    // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,AF,SO>::isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each row/column of the submatrix are guaranteed to
// conform to the alignment restrictions of the underlying element type.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool DenseSubmatrix<MT,AF,SO>::isAligned() const
{
   return isAligned_;
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
// rows and/or columns of the submatrix).
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool DenseSubmatrix<MT,AF,SO>::canSMPAssign() const
{
   return ( rows() > SMP_DMATASSIGN_THRESHOLD );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the column index (in case of
// a row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::IntrinsicType
   DenseSubmatrix<MT,AF,SO>::load( size_t i, size_t j ) const
{
   return loadu( i, j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the column index (in case of
// a row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DenseSubmatrix<MT,AF,SO>::IntrinsicType
   DenseSubmatrix<MT,AF,SO>::loadu( size_t i, size_t j ) const
{
   using blaze::load;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   if( isAligned_ ) {
      return matrix_.load( row_+i, column_+j );
   }
   else if( j != final_ ) {
      return matrix_.loadu( row_+i, column_+j );
   }
   else {
      AlignedArray<ElementType,IT::size> array;
      for( size_t k=0UL; k<rest_; ++k )
         array[k] = matrix_(row_+i,column_+j+k);
      for( size_t k=rest_; k<IT::size; ++k )
         array[k] = ElementType();
      return load( array.data() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store of a specific intrinsic element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the column index (in case of a row-major matrix) or the
// row index (in case of a column-major matrix) must be a multiple of the number of values inside
// the intrinsic element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void DenseSubmatrix<MT,AF,SO>::store( size_t i, size_t j, const IntrinsicType& value )
{
   storeu( i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the column index (in case of a row-major
// matrix) or the row index (in case of a column-major matrix) must be a multiple of the number
// of values inside the intrinsic element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void DenseSubmatrix<MT,AF,SO>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::store;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   if( isAligned_ ) {
      matrix_.store( row_+i, column_+j, value );
   }
   else if( j != final_ ) {
      matrix_.storeu( row_+i, column_+j, value );
   }
   else {
      AlignedArray<ElementType,IT::size> array;
      store( array.data(), value );
      for( size_t k=0UL; k<rest_; ++k )
         matrix_(row_+i,column_+j+k) = array[k];
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific intrinsic element of
// the dense submatrix. The row index must be smaller than the number of rows and the column
// index must be smaller than the number of columns. Additionally, the column index (in case
// of a row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void DenseSubmatrix<MT,AF,SO>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   storeu( i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,AF,SO>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,AF,SO>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row_+i,column_+j    ) = (~rhs)(i,j    );
         matrix_(row_+i,column_+j+1UL) = (~rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         matrix_(row_+i,column_+jpos) = (~rhs)(i,jpos);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,AF,SO>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,AF,SO>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   if( useStreaming && isAligned_ &&
       m_*n_ > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(~rhs).isAliased( &matrix_ ) )
   {
      for( size_t i=0UL; i<m_; ++i )
         for( size_t j=0UL; j<n_; j+=IT::size )
            matrix_.stream( row_+i, column_+j, (~rhs).load(i,j) );
   }
   else
   {
      const size_t jpos( n_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jpos, "Invalid end calculation" );

      for( size_t i=0UL; i<m_; ++i ) {
         typename MT2::ConstIterator it( (~rhs).begin(i) );
         for( size_t j=0UL; j<jpos; j+=IT::size*4UL ) {
            matrix_.storeu( row_+i, column_+j             , it.load() ); it += IT::size;
            matrix_.storeu( row_+i, column_+j+IT::size    , it.load() ); it += IT::size;
            matrix_.storeu( row_+i, column_+j+IT::size*2UL, it.load() ); it += IT::size;
            matrix_.storeu( row_+i, column_+j+IT::size*3UL, it.load() ); it += IT::size;
         }
         for( size_t j=jpos; j<n_; j+=IT::size, it+=IT::size ) {
            storeu( i, j, it.load() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,AF,SO>::assign( const DenseMatrix<MT2,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row_+i,column_+j) = (~rhs)(i,j);
            }
         }
      }
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
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,AF,SO>::assign( const SparseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) = element->value();
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
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,AF,SO>::assign( const SparseMatrix<MT2,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) = element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,AF,SO>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,AF,SO>::addAssign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row_+i,column_+j    ) += (~rhs)(i,j    );
         matrix_(row_+i,column_+j+1UL) += (~rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         matrix_(row_+i,column_+jpos) += (~rhs)(i,jpos);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,AF,SO>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,AF,SO>::addAssign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      typename MT2::ConstIterator it( (~rhs).begin(i) );
      for( size_t j=0UL; j<jpos; j+=IT::size*4UL ) {
         matrix_.storeu( row_+i, column_+j             , load(i,j             ) + it.load() ); it += IT::size;
         matrix_.storeu( row_+i, column_+j+IT::size    , load(i,j+IT::size    ) + it.load() ); it += IT::size;
         matrix_.storeu( row_+i, column_+j+IT::size*2UL, load(i,j+IT::size*2UL) + it.load() ); it += IT::size;
         matrix_.storeu( row_+i, column_+j+IT::size*3UL, load(i,j+IT::size*3UL) + it.load() ); it += IT::size;
      }
      for( size_t j=jpos; j<n_; j+=IT::size, it+=IT::size ) {
         storeu( i, j, load(i,j) + it.load() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,AF,SO>::addAssign( const DenseMatrix<MT2,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row_+i,column_+j) += (~rhs)(i,j);
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,AF,SO>::addAssign( const SparseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,AF,SO>::addAssign( const SparseMatrix<MT2,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,AF,SO>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,AF,SO>::subAssign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row_+i,column_+j    ) -= (~rhs)(i,j    );
         matrix_(row_+i,column_+j+1UL) -= (~rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         matrix_(row_+i,column_+jpos) -= (~rhs)(i,jpos);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,AF,SO>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,AF,SO>::subAssign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      typename MT2::ConstIterator it( (~rhs).begin(i) );
      for( size_t j=0UL; j<jpos; j+=IT::size*4UL ) {
         matrix_.storeu( row_+i, column_+j             , load(i,j             ) - it.load() ); it += IT::size;
         matrix_.storeu( row_+i, column_+j+IT::size    , load(i,j+IT::size    ) - it.load() ); it += IT::size;
         matrix_.storeu( row_+i, column_+j+IT::size*2UL, load(i,j+IT::size*2UL) - it.load() ); it += IT::size;
         matrix_.storeu( row_+i, column_+j+IT::size*3UL, load(i,j+IT::size*3UL) - it.load() ); it += IT::size;
      }
      for( size_t j=jpos; j<n_; j+=IT::size, it+=IT::size ) {
         storeu( i, j, load(i,j) - it.load() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,AF,SO>::subAssign( const DenseMatrix<MT2,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row_+i,column_+j) -= (~rhs)(i,j);
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,AF,SO>::subAssign( const SparseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) -= element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the dense matrix
        , bool AF         // Alignment flag
        , bool SO >       // Storage order
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,AF,SO>::subAssign( const SparseMatrix<MT2,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) -= element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR UNALIGNED COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubmatrix for unaligned column-major matrices.
// \ingroup dense_submatrix
//
// This specialization of DenseSubmatrix adapts the class template to the requirements of
// unaligned column-major matrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseSubmatrix<MT,unaligned,true> : public DenseMatrix< DenseSubmatrix<MT,unaligned,true>, true >
                                        , private Submatrix
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the matrix element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense submatrix will
       return references and iterators to const. Otherwise \a useConst will be set to 0 and
       the dense submatrix will offer write access to the dense matrix elements both via
       the function call operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubmatrix<MT,unaligned,true>   This;           //!< Type of this DenseSubmatrix instance.
   typedef typename SubmatrixTrait<MT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the submatrix elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the submatrix elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubmatrix&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant submatrix value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstPointer, ElementType* >::Type  Pointer;
   //**********************************************************************************************

   //**SubmatrixIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse submatrix.
   */
   template< typename IteratorType >  // Type of the dense matrix iterator
   class SubmatrixIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! The iterator category.
      typedef typename std::iterator_traits<IteratorType>::iterator_category  IteratorCategory;

      //! Type of the underlying elements.
      typedef typename std::iterator_traits<IteratorType>::value_type  ValueType;

      //! Pointer return type.
      typedef typename std::iterator_traits<IteratorType>::pointer  PointerType;

      //! Reference return type.
      typedef typename std::iterator_traits<IteratorType>::reference  ReferenceType;

      //! Difference between two iterators.
      typedef typename std::iterator_traits<IteratorType>::difference_type  DifferenceType;

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the SubmatrixIterator class.
      */
      inline SubmatrixIterator()
         : iterator_ (       )  // Iterator to the current submatrix element
         , final_    ( 0UL   )  // The final iterator for intrinsic operations
         , rest_     ( 0UL   )  // The number of remaining elements beyond the final iterator
         , isAligned_( false )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the SubmatrixIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param finalIterator The final iterator for intrinsic operations.
      // \param remainingElements The number of remaining elements beyond the final iterator.
      // \param isMemoryAligned Memory alignment flag.
      */
      inline SubmatrixIterator( IteratorType iterator, IteratorType finalIterator,
                                size_t remainingElements, bool isMemoryAligned )
         : iterator_ ( iterator          )  // Iterator to the current submatrix element
         , final_    ( finalIterator     )  // The final iterator for intrinsic operations
         , rest_     ( remainingElements )  // The number of remaining elements beyond the final iterator
         , isAligned_( isMemoryAligned   )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename IteratorType2 >
      inline SubmatrixIterator( const SubmatrixIterator<IteratorType2>& it )
         : iterator_ ( it.base()      )  // Iterator to the current submatrix element
         , final_    ( it.final()     )  // The final iterator for intrinsic operations
         , rest_     ( it.rest()      )  // The number of remaining elements beyond the final iterator
         , isAligned_( it.isAligned() )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline SubmatrixIterator& operator+=( size_t inc ) {
         iterator_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline SubmatrixIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubmatrixIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator++( int ) {
         return SubmatrixIterator( iterator_++, final_, rest_, isAligned_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline SubmatrixIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator--( int ) {
         return SubmatrixIterator( iterator_--, final_, rest_, isAligned_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReferenceType operator*() const {
         return *iterator_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Aligned load of an intrinsic element of the dense submatrix.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an aligned load of the current intrinsic element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType load() const {
         return loadu();
      }
      //*******************************************************************************************

      //**Loadu function***************************************************************************
      /*!\brief Unaligned load of an intrinsic element of the dense submatrix.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an unaligned load of the current intrinsic element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType loadu() const {
         if( isAligned_ ) {
            return iterator_.load();
         }
         else if( iterator_ != final_ ) {
            return iterator_.loadu();
         }
         else {
            AlignedArray<ElementType,IT::size> array;
            for( size_t i=0UL; i<rest_; ++i )
               array[i] = *(iterator_+i);
            for( size_t i=rest_; i<IT::size; ++i )
               array[i] = ElementType();
            return blaze::load( array.data() );
         }
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const SubmatrixIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const SubmatrixIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const SubmatrixIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const SubmatrixIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const SubmatrixIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const SubmatrixIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const SubmatrixIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( const SubmatrixIterator& it, size_t inc ) {
         return SubmatrixIterator( it.iterator_ + inc, it.final_, it.rest_, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a SubmatrixIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( size_t inc, const SubmatrixIterator& it ) {
         return SubmatrixIterator( it.iterator_ + inc, it.final_, it.rest_, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const SubmatrixIterator operator-( const SubmatrixIterator& it, size_t dec ) {
         return SubmatrixIterator( it.iterator_ - dec, it.final_, it.rest_, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the submatrix iterator.
      //
      // \return The current position of the submatrix iterator.
      */
      inline IteratorType base() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**Final function***************************************************************************
      /*!\brief Access to the final position of the submatrix iterator.
      //
      // \return The final position of the submatrix iterator.
      */
      inline IteratorType final() const {
         return final_;
      }
      //*******************************************************************************************

      //**Rest function****************************************************************************
      /*!\brief Access to the number of remaining elements beyond the final iterator.
      //
      // \return The number of remaining elements beyond the final iterator.
      */
      inline size_t rest() const {
         return rest_;
      }
      //*******************************************************************************************

      //**IsAligned function***********************************************************************
      /*!\brief Access to the iterator's memory alignment flag.
      //
      // \return \a true in case the iterator is aligned, \a false if it is not.
      */
      inline bool isAligned() const {
         return isAligned_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;   //!< Iterator to the current submatrix element.
      IteratorType final_;      //!< The final iterator for intrinsic operations.
      size_t       rest_;       //!< The number of remaining elements beyond the final iterator.
      bool         isAligned_;  //!< Memory alignment flag.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubmatrixIterator<typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename IfTrue< useConst, ConstIterator, SubmatrixIterator<typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = MT::vectorizable };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseSubmatrix& operator=( const ElementType& rhs );
   inline DenseSubmatrix& operator=( const DenseSubmatrix& rhs );

   template< typename MT2, bool SO >
   inline DenseSubmatrix& operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline DenseSubmatrix& operator*=( const Matrix<MT2,SO>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          rows() const;
                              inline size_t          columns() const;
                              inline size_t          spacing() const;
                              inline size_t          capacity() const;
                              inline size_t          capacity( size_t i ) const;
                              inline size_t          nonZeros() const;
                              inline size_t          nonZeros( size_t i ) const;
                              inline void            reset();
                              inline void            reset( size_t i );
                              inline DenseSubmatrix& transpose();
   template< typename Other > inline DenseSubmatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load ( size_t i, size_t j ) const;
   inline IntrinsicType loadu( size_t i, size_t j ) const;

   inline void store ( size_t i, size_t j, const IntrinsicType& value );
   inline void storeu( size_t i, size_t j, const IntrinsicType& value );
   inline void stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const;

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;   //!< The dense matrix containing the submatrix.
   const size_t row_;      //!< The first row of the submatrix.
   const size_t column_;   //!< The first column of the submatrix.
   const size_t m_;        //!< The number of rows of the submatrix.
   const size_t n_;        //!< The number of columns of the submatrix.
   const size_t rest_;     //!< The number of remaining elements in an unaligned intrinsic operation.
   const size_t final_;    //!< The final index for unaligned intrinsic operations.
                           /*!< In case the submatrix is not fully aligned and the submatrix is
                                involved in a vectorized operation, the final index indicates at
                                which index a special treatment for the remaining elements is
                                required. */
   const bool isAligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the submatrix is fully aligned.
                                In case the submatrix is fully aligned, no special handling has to
                                be used for the last elements of the submatrix in a vectorized
                                operation. In order to be aligned, the following conditions must
                                hold for the submatrix:
                                 - The first element of each row/column must be aligned
                                 - The submatrix must be at the end of the given matrix or
                                 - The number of rows/columns of the submatrix must be a multiple
                                   of the number of values per intrinsic element. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool AF2, bool SO2 > friend class DenseSubmatrix;

   template< bool AF1, typename MT2, bool AF2, bool SO2 >
   friend const DenseSubmatrix<MT2,AF1,SO2>
      submatrix( const DenseSubmatrix<MT2,AF2,SO2>& dm, size_t row, size_t column, size_t m, size_t n );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSymmetric( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isLower( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isUpper( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseMatrix<MT2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseMatrix<MT2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend typename DerestrictTrait< DenseSubmatrix<MT2,AF2,SO2> >::Type
      derestrict( DenseSubmatrix<MT2,AF2,SO2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
   BLAZE_STATIC_ASSERT( !IsRestricted<MT>::value || IsLower<MT>::value || IsUpper<MT>::value );
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
/*!\brief The constructor for DenseSubmatrix.
//
// \param matrix The dense matrix containing the submatrix.
// \param row The index of the first row of the submatrix in the given dense matrix.
// \param column The index of the first column of the submatrix in the given dense matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given dense matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,unaligned,true>::DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_   ( matrix       )  // The dense matrix containing the submatrix
   , row_      ( row          )  // The first row of the submatrix
   , column_   ( column       )  // The first column of the submatrix
   , m_        ( m            )  // The number of rows of the submatrix
   , n_        ( n            )  // The number of columns of the submatrix
   , rest_     ( m % IT::size )  // The number of remaining elements in an unaligned intrinsic operation
   , final_    ( m - rest_    )  // The final index for unaligned intrinsic operations
   , isAligned_( ( row % IT::size == 0UL ) &&
                 ( row + m == matrix.rows() || m % IT::size == 0UL ) )
{
   if( ( row + m > matrix.rows() ) || ( column + n > matrix.columns() ) )
      throw std::invalid_argument( "Invalid submatrix specification" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::Reference
   DenseSubmatrix<MT,unaligned,true>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::ConstReference
   DenseSubmatrix<MT,unaligned,true>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::Pointer DenseSubmatrix<MT,unaligned,true>::data()
{
   return matrix_.data() + row_ + column_*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::ConstPointer
   DenseSubmatrix<MT,unaligned,true>::data() const
{
   return matrix_.data() + row_ + column_*spacing();
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::Iterator
   DenseSubmatrix<MT,unaligned,true>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   const typename MT::Iterator first( matrix_.begin( column_ + j ) + row_ );
   return Iterator( first, first + final_, rest_, isAligned_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::ConstIterator
   DenseSubmatrix<MT,unaligned,true>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   const typename MT::ConstIterator first( matrix_.cbegin( column_ + j ) + row_ );
   return ConstIterator( first, first + final_, rest_, isAligned_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::ConstIterator
   DenseSubmatrix<MT,unaligned,true>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   const typename MT::ConstIterator first( matrix_.cbegin( column_ + j ) + row_ );
   return ConstIterator( first, first + final_, rest_, isAligned_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::Iterator
   DenseSubmatrix<MT,unaligned,true>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   const typename MT::Iterator last( matrix_.begin( column_ + j ) + row_ + m_ );
   return Iterator( last, last, rest_, isAligned_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::ConstIterator
   DenseSubmatrix<MT,unaligned,true>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   const typename MT::ConstIterator last( matrix_.cbegin( column_ + j ) + row_ + m_ );
   return ConstIterator( last, last, rest_, isAligned_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::ConstIterator
   DenseSubmatrix<MT,unaligned,true>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   const typename MT::ConstIterator last( matrix_.cbegin( column_ + j ) + row_ + m_ );
   return ConstIterator( last, last, rest_, isAligned_ );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,unaligned,true>&
   DenseSubmatrix<MT,unaligned,true>::operator=( const ElementType& rhs )
{
   const size_t jend( column_ + n_ );

   for( size_t j=column_; j<jend; ++j )
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

      for( size_t i=ibegin; i<iend; ++i )
         matrix_(i,j) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseSubmatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,unaligned,true>&
   DenseSubmatrix<MT,unaligned,true>::operator=( const DenseSubmatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ && column_ == rhs.column_ ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() )
      throw std::invalid_argument( "Submatrix sizes do not match" );

   if( !preservesInvariant( matrix_, rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

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
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,unaligned,true>&
   DenseSubmatrix<MT,unaligned,true>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   typedef typename If< IsAdaptor<MT>, typename MT2::CompositeType, const MT2& >::Type  Right;
   Right right( ~rhs );

   if( !preservesInvariant( matrix_, right ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   if( IsSparseMatrix<MT2>::value )
      reset();

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO  >     // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,unaligned,true>& >::Type
   DenseSubmatrix<MT,unaligned,true>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const AddType tmp( *this + (~rhs) );
      smpAssign( left, tmp );
   }
   else {
      smpAddAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO  >     // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,unaligned,true>& >::Type
   DenseSubmatrix<MT,unaligned,true>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const AddType tmp( *this + (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,unaligned,true>& >::Type
   DenseSubmatrix<MT,unaligned,true>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const SubType tmp( *this - (~rhs ) );
      smpAssign( left, tmp );
   }
   else {
      smpSubAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,unaligned,true>& >::Type
   DenseSubmatrix<MT,unaligned,true>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const SubType tmp( *this - (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,unaligned,true>&
   DenseSubmatrix<MT,unaligned,true>::operator*=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const MultType tmp( *this * (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a dense submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,unaligned,true> >::Type&
   DenseSubmatrix<MT,unaligned,true>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) * rhs );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense submatrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,unaligned,true> >::Type&
   DenseSubmatrix<MT,unaligned,true>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) / rhs );

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
/*!\brief Returns the number of rows of the dense submatrix.
//
// \return The number of rows of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::rows() const
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the dense submatrix.
//
// \return The number of columns of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::columns() const
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total
// number of elements of a column.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::spacing() const
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::capacity() const
{
   return rows() * columns();
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::capacity( size_t j ) const
{
   UNUSED_PARAMETER( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::nonZeros() const
{
   const size_t iend( row_ + m_ );
   const size_t jend( column_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t j=column_; j<jend; ++j )
      for( size_t i=row_; i<iend; ++i )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

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
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,unaligned,true>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( row_ + m_ );
   size_t nonzeros( 0UL );

   for( size_t i=row_; i<iend; ++i )
      if( !isDefault( matrix_(i,column_+j) ) )
         ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::reset()
{
   using blaze::clear;

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

      for( size_t i=ibegin; i<iend; ++i )
         clear( matrix_(i,j) );
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
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::reset( size_t j )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

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

   for( size_t i=ibegin; i<iend; ++i )
      clear( matrix_(i,column_+j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Transposing the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::runtime_error Invalid transpose of a non-quadratic submatrix.
// \exception std::runtime_error Invalid transpose of a lower matrix.
// \exception std::runtime_error Invalid transpose of an upper matrix.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// in case the underlying matrix is a lower or upper triangular matrix the function can only be
// used in case the submatrix does not contain elements from the upper or lower part of the matrix,
// respectively. The attempt to transpose a non-quadratic submatrix or an invalid part of a lower
// or triangular matrix results in a \a std::runtime_error exception.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,unaligned,true>& DenseSubmatrix<MT,unaligned,true>::transpose()
{
   if( rows() != columns() )
      throw std::runtime_error( "Invalid transpose of a non-quadratic submatrix" );

   if( IsLower<MT>::value && ( row_ + 1UL < column_ + n_ ) )
      throw std::runtime_error( "Invalid transpose of a lower matrix" );

   if( IsUpper<MT>::value && ( column_ + 1UL < row_ + m_ ) )
      throw std::runtime_error( "Invalid transpose of an upper matrix" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( trans(*this) );
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales all elements of the submatrix by the given scalar value \a scalar. Note
// that the function cannot be used to scale a submatrix on a lower or upper unitriangular matrix.
// The attempt to scale such a submatrix results in a compile time error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseSubmatrix<MT,unaligned,true>& DenseSubmatrix<MT,unaligned,true>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jend( column_ + n_ );

   for( size_t j=column_; j<jend; ++j )
   {
      const size_t ibegin( ( IsLower<MT>::value )
                           ?( ( IsStrictlyLower<MT>::value )
                              ?( max( j+1UL, row_ ) )
                              :( max( j, row_ ) ) )
                           :( row_ ) );
      const size_t iend  ( ( IsUpper<MT>::value )
                           ?( ( IsStrictlyUpper<MT>::value )
                              ?( min( j, row_+m_ ) )
                              :( min( j+1UL, row_+m_ ) ) )
                           :( row_+m_ ) );

      for( size_t i=ibegin; i<iend; ++i )
         matrix_(i,j) *= scalar;
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
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,unaligned,true>::hasOverlap() const
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric<MT>::value, "Unsymmetric matrix detected" );

   if( ( row_ + m_ <= column_ ) || ( column_ + n_ <= row_ ) )
      return false;
   else return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying matrix.
//
// \param lhs The matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying matrix of type \a MT would be
// violated by an assignment of the given matrix \a rhs. In case the matrix would be preserved,
// the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the left-hand side dense matrix
        , bool SO2       // Storage order of the left-hand side dense matrix
        , typename MT3   // Type of the right-hand side matrix
        , bool SO3 >     // Storage order of the right-hand side matrix
inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   UNUSED_PARAMETER( lhs, rhs );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying symmetric matrix.
//
// \param lhs The symmetric matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying symmetric matrix of type \a MT would
// be violated by an assignment of the given row-major matrix \a rhs. In case the matrix would be
// preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the left-hand side dense matrix
        , bool SO2       // Storage order of the left-hand side dense matrix
        , typename MT3   // Type of the right-hand side matrix
        , bool SO3 >     // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( !hasOverlap() )
      return true;

   const bool   lower( row_ > column_ );
   const size_t size ( min( row_ + m_, column_ + n_ ) - ( lower ? row_ : column_ ) );

   if( size < 2UL )
      return true;

   const size_t row   ( lower ? 0UL : column_ - row_ );
   const size_t column( lower ? row_ - column_ : 0UL );

   return isSymmetric( submatrix( ~rhs, row, column, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;

      const size_t jbegin( ( containsDiagonal )
                           ?( ( IsStrictlyLower<MT2>::value )
                              ?( row_ + i - column_ )
                              :( row_ + i - column_ + 1UL ) )
                           :( 0UL ) );

      for( size_t j=jbegin; j<n_; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t ipos( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                         ?( column_ + j - row_ + 1UL )
                         :( column_ + j - row_ ) );
      const size_t iend( min( ipos, m_ ) );

      for( size_t i=0UL; i<iend; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                             ?( row_ + i - column_ )
                             :( row_ + i - column_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(i) );
      RhsIterator element( (~rhs).lowerBound( i, index ) );

      if( IsUniLower<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != row_+i-column_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t index( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                          ?( column_ + j - row_ + 1UL )
                          :( column_ + j - row_ ) );
      const RhsIterator last( (~rhs).lowerBound( min( index, m_ ), j ) );

      if( IsUniLower<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(j) || ( last->index() != column_+j-row_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(j); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t jpos( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                         ?( row_ + i - column_ + 1UL )
                         :( row_ + i - column_ ) );
      const size_t jend( min( jpos, n_ ) );

      for( size_t j=0UL; j<jend; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;

      const size_t ibegin( ( containsDiagonal )
                           ?( ( IsStrictlyUpper<MT2>::value )
                              ?( column_ + j - row_ )
                              :( column_ + j - row_ + 1UL ) )
                           :( 0UL ) );

      for( size_t i=ibegin; i<m_; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t index( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                          ?( row_ + i - column_ + 1UL )
                          :( row_ + i - column_ ) );
      const RhsIterator last( (~rhs).lowerBound( i, min( index, n_ ) ) );

      if( IsUniUpper<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(i) || ( last->index() != row_+i-column_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(i); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                             ?( column_ + j - row_ )
                             :( column_ + j - row_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(j) );
      RhsIterator element( (~rhs).lowerBound( index, j ) );

      if( IsUniUpper<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != column_+j-row_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<n_; ++j ) {
         if( ( row_ + i != column_ + j ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<m_; ++i ) {
         if( ( column_ + j != row_ + i ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t i=0UL; i<m_; ++i ) {
      for( RhsIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element ) {
         if( ( row_ + i != column_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,unaligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t j=0UL; j<n_; ++j ) {
      for( RhsIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element ) {
         if( ( column_ + j != row_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
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
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,unaligned,true>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense submatrix
        , bool AF2       // Alignment flag of the foreign dense submatrix
        , bool SO2 >     // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,unaligned,true>::canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
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
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,unaligned,true>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense submatrix
        , bool AF2       // Alignment flag of the foreign dense submatrix
        , bool SO2 >     // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,unaligned,true>::isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each column of the submatrix are guaranteed to
// conform to the alignment restrictions of the underlying element type.
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,unaligned,true>::isAligned() const
{
   return isAligned_;
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
// rows and/or columns of the submatrix).
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,unaligned,true>::canSMPAssign() const
{
   return ( columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the intrinsic element. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::IntrinsicType
   DenseSubmatrix<MT,unaligned,true>::load( size_t i, size_t j ) const
{
   return loadu( i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the intrinsic element. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,unaligned,true>::IntrinsicType
   DenseSubmatrix<MT,unaligned,true>::loadu( size_t i, size_t j ) const
{
   using blaze::load;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   if( isAligned_ ) {
      return matrix_.load( row_+i, column_+j );
   }
   else if( i != final_ ) {
      return matrix_.loadu( row_+i, column_+j );
   }
   else {
      AlignedArray<ElementType,IT::size> array;
      for( size_t k=0UL; k<rest_; ++k )
         array[k] = matrix_(row_+i+k,column_+j);
      for( size_t k=rest_; k<IT::size; ++k )
         array[k] = ElementType();
      return load( array.data() );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store of a specific intrinsic element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of values
// inside the intrinsic element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::store( size_t i, size_t j, const IntrinsicType& value )
{
   storeu( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the row index must be a multiple of
// the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::store;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   if( isAligned_ ) {
      matrix_.store( row_+i, column_+j, value );
   }
   else if( i != final_ ) {
      matrix_.storeu( row_+i, column_+j, value );
   }
   else {
      AlignedArray<ElementType,IT::size> array;
      store( array.data(), value );
      for( size_t k=0UL; k<rest_; ++k )
         matrix_(row_+i+k,column_+j) = array[k];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific intrinsic element of
// the dense submatrix. The row index must be smaller than the number of rows and the column
// index must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the intrinsic element. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   storeu( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,unaligned,true>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,unaligned,true>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row_+i    ,column_+j) = (~rhs)(i    ,j);
         matrix_(row_+i+1UL,column_+j) = (~rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         matrix_(row_+ipos,column_+j) = (~rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,unaligned,true>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,unaligned,true>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   if( useStreaming && isAligned_ &&
       m_*n_ > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(~rhs).isAliased( &matrix_ ) )
   {
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=0UL; i<m_; i+=IT::size )
            matrix_.stream( row_+i, column_+j, (~rhs).load(i,j) );
   }
   else
   {
      const size_t ipos( m_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == ipos, "Invalid end calculation" );

      for( size_t j=0UL; j<n_; ++j ) {
         typename MT2::ConstIterator it( (~rhs).begin(j) );
         for( size_t i=0UL; i<ipos; i+=IT::size*4UL ) {
            matrix_.storeu( row_+i             , column_+j, it.load() ); it += IT::size;
            matrix_.storeu( row_+i+IT::size    , column_+j, it.load() ); it += IT::size;
            matrix_.storeu( row_+i+IT::size*2UL, column_+j, it.load() ); it += IT::size;
            matrix_.storeu( row_+i+IT::size*3UL, column_+j, it.load() ); it += IT::size;
         }
         for( size_t i=ipos; i<m_; i+=IT::size, it+=IT::size ) {
            storeu( i, j, it.load() );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row_+i,column_+j) = (~rhs)(i,j);
            }
         }
      }
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
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,unaligned,true>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) = element->value();
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
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,unaligned,true>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,unaligned,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,unaligned,true>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row_+i    ,column_+j) += (~rhs)(i    ,j);
         matrix_(row_+i+1UL,column_+j) += (~rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         matrix_(row_+ipos,column_+j) += (~rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,unaligned,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,unaligned,true>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      typename MT2::ConstIterator it( (~rhs).begin(j) );
      for( size_t i=0UL; i<ipos; i+=IT::size*4UL ) {
         matrix_.storeu( row_+i             , column_+j, load(i             ,j) + it.load() ); it += IT::size;
         matrix_.storeu( row_+i+IT::size    , column_+j, load(i+IT::size    ,j) + it.load() ); it += IT::size;
         matrix_.storeu( row_+i+IT::size*2UL, column_+j, load(i+IT::size*2UL,j) + it.load() ); it += IT::size;
         matrix_.storeu( row_+i+IT::size*3UL, column_+j, load(i+IT::size*3UL,j) + it.load() ); it += IT::size;
      }
      for( size_t i=ipos; i<m_; i+=IT::size, it+=IT::size ) {
         storeu( i, j, load(i,j) + it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row_+i,column_+j) += (~rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,unaligned,true>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,unaligned,true>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,unaligned,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,unaligned,true>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row_+i    ,column_+j) -= (~rhs)(i    ,j);
         matrix_(row_+i+1UL,column_+j) -= (~rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         matrix_(row_+ipos,column_+j) -= (~rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,unaligned,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,unaligned,true>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      typename MT2::ConstIterator it( (~rhs).begin(j) );
      for( size_t i=0UL; i<ipos; i+=IT::size*4UL ) {
         matrix_.storeu( row_+i             , column_+j, load(i             ,j) - it.load() ); it += IT::size;
         matrix_.storeu( row_+i+IT::size    , column_+j, load(i+IT::size    ,j) - it.load() ); it += IT::size;
         matrix_.storeu( row_+i+IT::size*2UL, column_+j, load(i+IT::size*2UL,j) - it.load() ); it += IT::size;
         matrix_.storeu( row_+i+IT::size*3UL, column_+j, load(i+IT::size*3UL,j) - it.load() ); it += IT::size;
      }
      for( size_t i=ipos; i<m_; i+=IT::size, it+=IT::size ) {
         storeu( i, j, load(i,j) - it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,unaligned,true>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row_+i,column_+j) -= (~rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,unaligned,true>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,unaligned,true>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ALIGNED ROW-MAJOR SUBMATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubmatrix for aligned row-major submatrices.
// \ingroup dense_submatrix
//
// This specialization of DenseSubmatrix adapts the class template to the requirements of
// aligned row-major submatrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseSubmatrix<MT,aligned,false> : public DenseMatrix< DenseSubmatrix<MT,aligned,false>, false >
                                       , private Submatrix
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the matrix element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense submatrix will
       return references and iterators to const. Otherwise \a useConst will be set to 0 and
       the dense submatrix will offer write access to the dense matrix elements both via
       the function call operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubmatrix<MT,aligned,false>    This;           //!< Type of this DenseSubmatrix instance.
   typedef typename SubmatrixTrait<MT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the submatrix elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the submatrix elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubmatrix&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant submatrix value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename IfTrue< useConst, ConstIterator, typename MT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = MT::vectorizable };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseSubmatrix& operator=( const ElementType& rhs );
   inline DenseSubmatrix& operator=( const DenseSubmatrix& rhs );

   template< typename MT2, bool SO >
   inline DenseSubmatrix& operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline DenseSubmatrix& operator*=( const Matrix<MT2,SO>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          rows() const;
                              inline size_t          columns() const;
                              inline size_t          spacing() const;
                              inline size_t          capacity() const;
                              inline size_t          capacity( size_t i ) const;
                              inline size_t          nonZeros() const;
                              inline size_t          nonZeros( size_t i ) const;
                              inline void            reset();
                              inline void            reset( size_t i );
                              inline DenseSubmatrix& transpose();
   template< typename Other > inline DenseSubmatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   BLAZE_ALWAYS_INLINE IntrinsicType load ( size_t i, size_t j ) const;
   BLAZE_ALWAYS_INLINE IntrinsicType loadu( size_t i, size_t j ) const;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,false>& rhs );

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const;

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The dense matrix containing the submatrix.
   const size_t row_;     //!< The first row of the submatrix.
   const size_t column_;  //!< The first column of the submatrix.
   const size_t m_;       //!< The number of rows of the submatrix.
   const size_t n_;       //!< The number of columns of the submatrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool AF2, bool SO2 > friend class DenseSubmatrix;

   template< bool AF1, typename MT2, bool AF2, bool SO2 >
   friend const DenseSubmatrix<MT2,AF1,SO2>
      submatrix( const DenseSubmatrix<MT2,AF2,SO2>& dm, size_t row, size_t column, size_t m, size_t n );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSymmetric( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isLower( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isUpper( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseMatrix<MT2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseMatrix<MT2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend typename DerestrictTrait< DenseSubmatrix<MT2,AF2,SO2> >::Type
      derestrict( DenseSubmatrix<MT2,AF2,SO2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
   BLAZE_STATIC_ASSERT( !IsRestricted<MT>::value || IsLower<MT>::value || IsUpper<MT>::value );
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
/*!\brief The constructor for DenseSubmatrix.
//
// \param matrix The dense matrix containing the submatrix.
// \param row The index of the first row of the submatrix in the given dense matrix.
// \param column The index of the first column of the submatrix in the given dense matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given dense matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,false>::DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_( matrix )  // The dense matrix containing the submatrix
   , row_   ( row    )  // The first row of the submatrix
   , column_( column )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
{
   if( ( row + m > matrix.rows() ) || ( column + n > matrix.columns() ) )
      throw std::invalid_argument( "Invalid submatrix specification" );

   if( column % IT::size != 0UL || ( column_ + n_ != matrix_.columns() && n_ % IT::size != 0UL ) )
      throw std::invalid_argument( "Invalid submatrix alignment" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::Reference
   DenseSubmatrix<MT,aligned,false>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::ConstReference
   DenseSubmatrix<MT,aligned,false>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::Pointer DenseSubmatrix<MT,aligned,false>::data()
{
   return matrix_.data() + row_*spacing() + column_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::ConstPointer
   DenseSubmatrix<MT,aligned,false>::data() const
{
   return matrix_.data() + row_*spacing() + column_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::Iterator
   DenseSubmatrix<MT,aligned,false>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.begin( row_ + i ) + column_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::ConstIterator
   DenseSubmatrix<MT,aligned,false>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row_ + i ) + column_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::ConstIterator
   DenseSubmatrix<MT,aligned,false>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row_ + i ) + column_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::Iterator
   DenseSubmatrix<MT,aligned,false>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.begin( row_ + i ) + column_ + n_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::ConstIterator
   DenseSubmatrix<MT,aligned,false>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row_ + i ) + column_ + n_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,false>::ConstIterator
   DenseSubmatrix<MT,aligned,false>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row_ + i ) + column_ + n_ );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,false>&
   DenseSubmatrix<MT,aligned,false>::operator=( const ElementType& rhs )
{
   const size_t iend( row_ + m_ );

   for( size_t i=row_; i<iend; ++i )
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

      for( size_t j=jbegin; j<jend; ++j )
         matrix_(i,j) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseSubmatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,false>&
   DenseSubmatrix<MT,aligned,false>::operator=( const DenseSubmatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ && column_ == rhs.column_ ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() )
      throw std::invalid_argument( "Submatrix sizes do not match" );

   if( !preservesInvariant( matrix_, rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

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
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,aligned,false>&
   DenseSubmatrix<MT,aligned,false>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   typedef typename If< IsAdaptor<MT>, typename MT2::CompositeType, const MT2& >::Type  Right;
   Right right( ~rhs );

   if( !preservesInvariant( matrix_, right ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   if( IsSparseMatrix<MT2>::value )
      reset();

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,aligned,false>& >::Type
   DenseSubmatrix<MT,aligned,false>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const AddType tmp( *this + (~rhs) );
      smpAssign( left, tmp );
   }
   else {
      smpAddAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,aligned,false>& >::Type
   DenseSubmatrix<MT,aligned,false>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const AddType tmp( *this + (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,aligned,false>& >::Type
   DenseSubmatrix<MT,aligned,false>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const SubType tmp( *this - (~rhs ) );
      smpAssign( left, tmp );
   }
   else {
      smpSubAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,aligned,false>& >::Type
   DenseSubmatrix<MT,aligned,false>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const SubType tmp( *this - (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,aligned,false>&
   DenseSubmatrix<MT,aligned,false>::operator*=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const MultType tmp( *this * (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a dense submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,aligned,false> >::Type&
   DenseSubmatrix<MT,aligned,false>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) * rhs );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense submatrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,aligned,false> >::Type&
   DenseSubmatrix<MT,aligned,false>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) / rhs );

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
/*!\brief Returns the number of rows of the dense submatrix.
//
// \return The number of rows of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::rows() const
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the dense submatrix.
//
// \return The number of columns of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::columns() const
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two rows/columns.
//
// \return The spacing between the beginning of two rows/columns.
//
// This function returns the spacing between the beginning of two rows/columns, i.e. the
// total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::spacing() const
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::capacity() const
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::capacity( size_t i ) const
{
   UNUSED_PARAMETER( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::nonZeros() const
{
   const size_t iend( row_ + m_ );
   const size_t jend( column_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t i=row_; i<iend; ++i )
      for( size_t j=column_; j<jend; ++j )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,false>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( column_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t j=column_; j<jend; ++j )
      if( !isDefault( matrix_(row_+i,j) ) )
         ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,aligned,false>::reset()
{
   using blaze::clear;

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

      for( size_t j=jbegin; j<jend; ++j )
         clear( matrix_(i,j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,aligned,false>::reset( size_t i )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

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

   for( size_t j=jbegin; j<jend; ++j )
      clear( matrix_(row_+i,j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Transposing the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::runtime_error Invalid transpose of a non-quadratic submatrix.
// \exception std::runtime_error Invalid transpose of a lower matrix.
// \exception std::runtime_error Invalid transpose of an upper matrix.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// in case the underlying matrix is a lower or upper triangular matrix the function can only be
// used in case the submatrix does not contain elements from the upper or lower part of the matrix,
// respectively. The attempt to transpose a non-quadratic submatrix or an invalid part of a lower
// or triangular matrix results in a \a std::runtime_error exception.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,false>& DenseSubmatrix<MT,aligned,false>::transpose()
{
   if( rows() != columns() )
      throw std::runtime_error( "Invalid transpose of a non-quadratic submatrix" );

   if( IsLower<MT>::value && ( row_ + 1UL < column_ + n_ ) )
      throw std::runtime_error( "Invalid transpose of a lower matrix" );

   if( IsUpper<MT>::value && ( column_ + 1UL < row_ + m_ ) )
      throw std::runtime_error( "Invalid transpose of an upper matrix" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( trans(*this) );
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales all elements of the submatrix by the given scalar value \a scalar. Note
// that the function cannot be used to scale a submatrix on a lower or upper unitriangular matrix.
// The attempt to scale such a submatrix results in a compile time error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseSubmatrix<MT,aligned,false>& DenseSubmatrix<MT,aligned,false>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t iend( row_ + m_ );

   for( size_t i=row_; i<iend; ++i )
   {
      const size_t jbegin( ( IsUpper<MT>::value )
                           ?( ( IsStrictlyUpper<MT>::value )
                              ?( max( i+1UL, column_ ) )
                              :( max( i, column_ ) ) )
                           :( column_ ) );
      const size_t jend  ( ( IsLower<MT>::value )
                           ?( ( IsStrictlyLower<MT>::value )
                              ?( min( i, column_+n_ ) )
                              :( min( i+1UL, column_+n_ ) ) )
                           :( column_+n_ ) );

      for( size_t j=jbegin; j<jend; ++j )
         matrix_(i,j) *= scalar;
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
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,aligned,false>::hasOverlap() const
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric<MT>::value, "Unsymmetric matrix detected" );

   if( ( row_ + m_ <= column_ ) || ( column_ + n_ <= row_ ) )
      return false;
   else return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying matrix.
//
// \param lhs The matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying matrix of type \a MT would be
// violated by an assignment of the given matrix \a rhs. In case the matrix would be preserved,
// the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the left-hand side dense matrix
        , bool SO2       // Storage order of the left-hand side dense matrix
        , typename MT3   // Type of the right-hand side matrix
        , bool SO3 >     // Storage order of the right-hand side matrix
inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   UNUSED_PARAMETER( lhs, rhs );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying symmetric matrix.
//
// \param lhs The symmetric matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying symmetric matrix of type \a MT would
// be violated by an assignment of the given row-major matrix \a rhs. In case the matrix would be
// preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the left-hand side dense matrix
        , bool SO2       // Storage order of the left-hand side dense matrix
        , typename MT3   // Type of the right-hand side matrix
        , bool SO3 >     // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( !hasOverlap() )
      return true;

   const bool   lower( row_ > column_ );
   const size_t size ( min( row_ + m_, column_ + n_ ) - ( lower ? row_ : column_ ) );

   if( size < 2UL )
      return true;

   const size_t row   ( lower ? 0UL : column_ - row_ );
   const size_t column( lower ? row_ - column_ : 0UL );

   return isSymmetric( submatrix( ~rhs, row, column, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;

      const size_t jbegin( ( containsDiagonal )
                           ?( ( IsStrictlyLower<MT2>::value )
                              ?( row_ + i - column_ )
                              :( row_ + i - column_ + 1UL ) )
                           :( 0UL ) );

      for( size_t j=jbegin; j<n_; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t ipos( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                         ?( column_ + j - row_ + 1UL )
                         :( column_ + j - row_ ) );
      const size_t iend( min( ipos, m_ ) );

      for( size_t i=0UL; i<iend; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                             ?( row_ + i - column_ )
                             :( row_ + i - column_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(i) );
      RhsIterator element( (~rhs).lowerBound( i, index ) );

      if( IsUniLower<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != row_+i-column_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t index( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                          ?( column_ + j - row_ + 1UL )
                          :( column_ + j - row_ ) );
      const RhsIterator last( (~rhs).lowerBound( min( index, m_ ), j ) );

      if( IsUniLower<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(j) || ( last->index() != column_+j-row_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(j); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t jpos( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                         ?( row_ + i - column_ + 1UL )
                         :( row_ + i - column_ ) );
      const size_t jend( min( jpos, n_ ) );

      for( size_t j=0UL; j<jend; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;

      const size_t ibegin( ( containsDiagonal )
                           ?( ( IsStrictlyUpper<MT2>::value )
                              ?( column_ + j - row_ )
                              :( column_ + j - row_ + 1UL ) )
                           :( 0UL ) );

      for( size_t i=ibegin; i<m_; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t index( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                          ?( row_ + i - column_ + 1UL )
                          :( row_ + i - column_ ) );
      const RhsIterator last( (~rhs).lowerBound( i, min( index, n_ ) ) );

      if( IsUniUpper<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(i) || ( last->index() != row_+i-column_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(i); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                             ?( column_ + j - row_ )
                             :( column_ + j - row_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(j) );
      RhsIterator element( (~rhs).lowerBound( index, j ) );

      if( IsUniUpper<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != column_+j-row_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<n_; ++j ) {
         if( ( row_ + i != column_ + j ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<m_; ++i ) {
         if( ( column_ + j != row_ + i ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t i=0UL; i<m_; ++i ) {
      for( RhsIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element ) {
         if( ( row_ + i != column_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,false>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t j=0UL; j<n_; ++j ) {
      for( RhsIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element ) {
         if( ( column_ + j != row_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
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
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,aligned,false>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense submatrix
        , bool AF2       // Alignment flag of the foreign dense submatrix
        , bool SO2 >     // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,aligned,false>::canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
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
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,aligned,false>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense submatrix
        , bool AF2       // Alignment flag of the foreign dense submatrix
        , bool SO2 >     // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,aligned,false>::isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each row of the submatrix are guaranteed to conform
// to the alignment restrictions of the underlying element type.
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,aligned,false>::isAligned() const
{
   return true;
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
// rows and/or columns of the submatrix).
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,aligned,false>::canSMPAssign() const
{
   return ( rows() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the column index (in case of
// a row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseSubmatrix<MT,aligned,false>::IntrinsicType
   DenseSubmatrix<MT,aligned,false>::load( size_t i, size_t j ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   return matrix_.load( row_+i, column_+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the column index (in case of
// a row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseSubmatrix<MT,aligned,false>::IntrinsicType
   DenseSubmatrix<MT,aligned,false>::loadu( size_t i, size_t j ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   return matrix_.loadu( row_+i, column_+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store of a specific intrinsic element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the column index (in case of a row-major matrix) or the
// row index (in case of a column-major matrix) must be a multiple of the number of values inside
// the intrinsic element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void
   DenseSubmatrix<MT,aligned,false>::store( size_t i, size_t j, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   return matrix_.store( row_+i, column_+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the column index (in case of a row-major
// matrix) or the row index (in case of a column-major matrix) must be a multiple of the number
// of values inside the intrinsic element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void
   DenseSubmatrix<MT,aligned,false>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   matrix_.storeu( row_+i, column_+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific intrinsic element of
// the dense submatrix. The row index must be smaller than the number of rows and the column
// index must be smaller than the number of columns. Additionally, the column index (in case
// of a row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void
   DenseSubmatrix<MT,aligned,false>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   matrix_.stream( row_+i, column_+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,aligned,false>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,false>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row_+i,column_+j    ) = (~rhs)(i,j    );
         matrix_(row_+i,column_+j+1UL) = (~rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         matrix_(row_+i,column_+jpos) = (~rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,aligned,false>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,false>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   if( useStreaming &&
       m_*n_ > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(~rhs).isAliased( &matrix_ ) )
   {
      for( size_t i=0UL; i<m_; ++i )
         for( size_t j=0UL; j<n_; j+=IT::size )
            stream( i, j, (~rhs).load(i,j) );
   }
   else
   {
      const size_t jpos( n_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jpos, "Invalid end calculation" );

      for( size_t i=0UL; i<m_; ++i ) {
         typename MT2::ConstIterator it( (~rhs).begin(i) );
         for( size_t j=0UL; j<jpos; j+=IT::size*4UL ) {
            store( i, j             , it.load() ); it += IT::size;
            store( i, j+IT::size    , it.load() ); it += IT::size;
            store( i, j+IT::size*2UL, it.load() ); it += IT::size;
            store( i, j+IT::size*3UL, it.load() ); it += IT::size;
         }
         for( size_t j=jpos; j<n_; j+=IT::size, it+=IT::size ) {
            store( i, j, it.load() );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,aligned,false>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row_+i,column_+j) = (~rhs)(i,j);
            }
         }
      }
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
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,false>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) = element->value();
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
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,false>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,aligned,false>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,false>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row_+i,column_+j    ) += (~rhs)(i,j    );
         matrix_(row_+i,column_+j+1UL) += (~rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         matrix_(row_+i,column_+jpos) += (~rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,aligned,false>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,false>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      typename MT2::ConstIterator it( (~rhs).begin(i) );
      for( size_t j=0UL; j<jpos; j+=IT::size*4UL ) {
         store( i, j             , load(i,j             ) + it.load() ); it += IT::size;
         store( i, j+IT::size    , load(i,j+IT::size    ) + it.load() ); it += IT::size;
         store( i, j+IT::size*2UL, load(i,j+IT::size*2UL) + it.load() ); it += IT::size;
         store( i, j+IT::size*3UL, load(i,j+IT::size*3UL) + it.load() ); it += IT::size;
      }
      for( size_t j=jpos; j<n_; j+=IT::size, it+=IT::size ) {
         store( i, j, load(i,j) + it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,aligned,false>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row_+i,column_+j) += (~rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,false>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,false>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,aligned,false>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,false>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row_+i,column_+j    ) -= (~rhs)(i,j    );
         matrix_(row_+i,column_+j+1UL) -= (~rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         matrix_(row_+i,column_+jpos) -= (~rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,aligned,false>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,false>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jpos( n_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jpos, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      typename MT2::ConstIterator it( (~rhs).begin(i) );
      for( size_t j=0UL; j<jpos; j+=IT::size*4UL ) {
         store( i, j             , load(i,j             ) - it.load() ); it += IT::size;
         store( i, j+IT::size    , load(i,j+IT::size    ) - it.load() ); it += IT::size;
         store( i, j+IT::size*2UL, load(i,j+IT::size*2UL) - it.load() ); it += IT::size;
         store( i, j+IT::size*3UL, load(i,j+IT::size*3UL) - it.load() ); it += IT::size;
      }
      for( size_t j=jpos; j<n_; j+=IT::size, it+=IT::size ) {
         store( i, j, load(i,j) - it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,aligned,false>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row_+i,column_+j) -= (~rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,false>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,false>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ALIGNED COLUMN-MAJOR SUBMATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubmatrix for aligned column-major matrices.
// \ingroup dense_submatrix
//
// This specialization of DenseSubmatrix adapts the class template to the requirements of
// aligned column-major matrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseSubmatrix<MT,aligned,true> : public DenseMatrix< DenseSubmatrix<MT,aligned,true>, true >
                                      , private Submatrix
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the matrix element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense submatrix will
       return references and iterators to const. Otherwise \a useConst will be set to 0 and
       the dense submatrix will offer write access to the dense matrix elements both via
       the function call operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubmatrix<MT,aligned,true>     This;           //!< Type of this DenseSubmatrix instance.
   typedef typename SubmatrixTrait<MT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the submatrix elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the submatrix elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubmatrix&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant submatrix value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant submatrix value.
   typedef typename IfTrue< useConst, ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename IfTrue< useConst, ConstIterator, typename MT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = MT::vectorizable };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseSubmatrix& operator=( const ElementType& rhs );
   inline DenseSubmatrix& operator=( const DenseSubmatrix& rhs );

   template< typename MT2, bool SO >
   inline DenseSubmatrix& operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator+=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >, DenseSubmatrix& >::Type
      operator-=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline DenseSubmatrix& operator*=( const Matrix<MT2,SO>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          rows() const;
                              inline size_t          columns() const;
                              inline size_t          spacing() const;
                              inline size_t          capacity() const;
                              inline size_t          capacity( size_t i ) const;
                              inline size_t          nonZeros() const;
                              inline size_t          nonZeros( size_t i ) const;
                              inline void            reset();
                              inline void            reset( size_t i );
                              inline DenseSubmatrix& transpose();
   template< typename Other > inline DenseSubmatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT2::vectorizable &&
                     IsSame<ElementType,typename MT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool AF2, bool SO2 >
   inline bool isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   BLAZE_ALWAYS_INLINE IntrinsicType load ( size_t i, size_t j ) const;
   BLAZE_ALWAYS_INLINE IntrinsicType loadu( size_t i, size_t j ) const;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAssign<MT2> >::Type
      assign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedAddAssign<MT2> >::Type
      addAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline typename DisableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline typename EnableIf< VectorizedSubAssign<MT2> >::Type
      subAssign( const DenseMatrix<MT2,true>& rhs );

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const;

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3, bool SO3 >
   inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs );

   template< typename MT2, bool SO2, typename MT3 >
   inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
      preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The dense matrix containing the submatrix.
   const size_t row_;     //!< The first row of the submatrix.
   const size_t column_;  //!< The first column of the submatrix.
   const size_t m_;       //!< The number of rows of the submatrix.
   const size_t n_;       //!< The number of columns of the submatrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool AF2, bool SO2 > friend class DenseSubmatrix;

   template< bool AF1, typename MT2, bool AF2, bool SO2 >
   friend const DenseSubmatrix<MT2,AF1,SO2>
      submatrix( const DenseSubmatrix<MT2,AF2,SO2>& dm, size_t row, size_t column, size_t m, size_t n );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSymmetric( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isLower( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isUpper( const DenseSubmatrix<MT2,AF2,SO2>& dm );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseMatrix<MT2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseMatrix<MT2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend bool isSame( const DenseSubmatrix<MT2,AF2,SO2>& a, const DenseSubmatrix<MT2,AF2,SO2>& b );

   template< typename MT2, bool AF2, bool SO2 >
   friend typename DerestrictTrait< DenseSubmatrix<MT2,AF2,SO2> >::Type
      derestrict( DenseSubmatrix<MT2,AF2,SO2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
   BLAZE_STATIC_ASSERT( !IsRestricted<MT>::value || IsLower<MT>::value || IsUpper<MT>::value );
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
/*!\brief The constructor for DenseSubmatrix.
//
// \param matrix The dense matrix containing the submatrix.
// \param row The index of the first row of the submatrix in the given dense matrix.
// \param column The index of the first column of the submatrix in the given dense matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given dense matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,true>::DenseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_( matrix )  // The dense matrix containing the submatrix
   , row_   ( row    )  // The first row of the submatrix
   , column_( column )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
{
   if( ( row + m > matrix.rows() ) || ( column + n > matrix.columns() ) )
      throw std::invalid_argument( "Invalid submatrix specification" );

   if( row % IT::size != 0UL || ( row_ + m_ != matrix_.rows() && m_ % IT::size != 0UL ) )
      throw std::invalid_argument( "Invalid submatrix alignment" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::Reference
   DenseSubmatrix<MT,aligned,true>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::ConstReference
   DenseSubmatrix<MT,aligned,true>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row_+i,column_+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::Pointer DenseSubmatrix<MT,aligned,true>::data()
{
   return matrix_.data() + row_ + column_*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::ConstPointer
   DenseSubmatrix<MT,aligned,true>::data() const
{
   return matrix_.data() + row_ + column_*spacing();
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::Iterator
   DenseSubmatrix<MT,aligned,true>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.begin( column_ + j ) + row_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::ConstIterator
   DenseSubmatrix<MT,aligned,true>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column_ + j ) + row_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::ConstIterator
   DenseSubmatrix<MT,aligned,true>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column_ + j ) + row_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::Iterator
   DenseSubmatrix<MT,aligned,true>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.begin( column_ + j ) + row_ + m_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::ConstIterator
   DenseSubmatrix<MT,aligned,true>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column_ + j ) + row_ + m_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseSubmatrix<MT,aligned,true>::ConstIterator
   DenseSubmatrix<MT,aligned,true>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column_ + j ) + row_ + m_ );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,true>&
   DenseSubmatrix<MT,aligned,true>::operator=( const ElementType& rhs )
{
   const size_t jend( column_ + n_ );

   for( size_t j=column_; j<jend; ++j )
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

      for( size_t i=ibegin; i<iend; ++i )
         matrix_(i,j) = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseSubmatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,true>&
   DenseSubmatrix<MT,aligned,true>::operator=( const DenseSubmatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row_ == rhs.row_ && column_ == rhs.column_ ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() )
      throw std::invalid_argument( "Submatrix sizes do not match" );

   if( !preservesInvariant( matrix_, rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

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
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the
// current sizes of the two matrices don't match, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is a symmetric matrix and the assignment
// would violate its symmetry, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,aligned,true>&
   DenseSubmatrix<MT,aligned,true>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   typedef typename If< IsAdaptor<MT>, typename MT2::CompositeType, const MT2& >::Type  Right;
   Right right( ~rhs );

   if( !preservesInvariant( matrix_, right ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   if( IsSparseMatrix<MT2>::value )
      reset();

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO  >     // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,aligned,true>& >::Type
   DenseSubmatrix<MT,aligned,true>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const AddType tmp( *this + (~rhs) );
      smpAssign( left, tmp );
   }
   else {
      smpAddAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO  >     // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,aligned,true>& >::Type
   DenseSubmatrix<MT,aligned,true>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename AddTrait<ResultType,typename MT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const AddType tmp( *this + (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename DisableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                         , DenseSubmatrix<MT,aligned,true>& >::Type
   DenseSubmatrix<MT,aligned,true>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( !preservesInvariant( matrix_, ~rhs ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( ( IsSymmetric<MT>::value && hasOverlap() ) || (~rhs).canAlias( &matrix_ ) ) {
      const SubType tmp( *this - (~rhs ) );
      smpAssign( left, tmp );
   }
   else {
      smpSubAssign( left, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsAdaptor<MT>, RequiresEvaluation<MT2> >
                        , DenseSubmatrix<MT,aligned,true>& >::Type
   DenseSubmatrix<MT,aligned,true>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename SubTrait<ResultType,typename MT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const SubType tmp( *this - (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to matrix adaptor.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline DenseSubmatrix<MT,aligned,true>&
   DenseSubmatrix<MT,aligned,true>::operator*=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const MultType tmp( *this * (~rhs) );

   if( !preservesInvariant( matrix_, tmp ) )
      throw std::invalid_argument( "Invalid assignment to matrix adaptor" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( !IsLower<MT>::value || isLower( derestrict( matrix_ ) ), "Lower violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsUpper<MT>::value || isUpper( derestrict( matrix_ ) ), "Upper violation detected" );
   BLAZE_INTERNAL_ASSERT( !IsSymmetric<MT>::value || isSymmetric( derestrict( matrix_ ) ), "Symmetry violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a dense submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,aligned,true> >::Type&
   DenseSubmatrix<MT,aligned,true>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) * rhs );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense submatrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the dense submatrix.
//
// This operator cannot be used for submatrices on lower or upper unitriangular matrices. The
// attempt to scale such a submatrix results in a compilation error!
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubmatrix<MT,aligned,true> >::Type&
   DenseSubmatrix<MT,aligned,true>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) / rhs );

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
/*!\brief Returns the number of rows of the dense submatrix.
//
// \return The number of rows of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::rows() const
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the dense submatrix.
//
// \return The number of columns of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::columns() const
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total
// number of elements of a column.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::spacing() const
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::capacity() const
{
   return rows() * columns();
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::capacity( size_t j ) const
{
   UNUSED_PARAMETER( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::nonZeros() const
{
   const size_t iend( row_ + m_ );
   const size_t jend( column_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t j=column_; j<jend; ++j )
      for( size_t i=row_; i<iend; ++i )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

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
template< typename MT >  // Type of the dense matrix
inline size_t DenseSubmatrix<MT,aligned,true>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( row_ + m_ );
   size_t nonzeros( 0UL );

   for( size_t i=row_; i<iend; ++i )
      if( !isDefault( matrix_(i,column_+j) ) )
         ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,aligned,true>::reset()
{
   using blaze::clear;

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

      for( size_t i=ibegin; i<iend; ++i )
         clear( matrix_(i,j) );
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
template< typename MT >  // Type of the dense matrix
inline void DenseSubmatrix<MT,aligned,true>::reset( size_t j )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

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

   for( size_t i=ibegin; i<iend; ++i )
      clear( matrix_(i,column_+j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Transposing the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::runtime_error Invalid transpose of a non-quadratic submatrix.
// \exception std::runtime_error Invalid transpose of a lower matrix.
// \exception std::runtime_error Invalid transpose of an upper matrix.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// in case the underlying matrix is a lower or upper triangular matrix the function can only be
// used in case the submatrix does not contain elements from the upper or lower part of the matrix,
// respectively. The attempt to transpose a non-quadratic submatrix or an invalid part of a lower
// or triangular matrix results in a \a std::runtime_error exception.
*/
template< typename MT >  // Type of the dense matrix
inline DenseSubmatrix<MT,aligned,true>& DenseSubmatrix<MT,aligned,true>::transpose()
{
   if( rows() != columns() )
      throw std::runtime_error( "Invalid transpose of a non-quadratic submatrix" );

   if( IsLower<MT>::value && ( row_ + 1UL < column_ + n_ ) )
      throw std::runtime_error( "Invalid transpose of a lower matrix" );

   if( IsUpper<MT>::value && ( column_ + 1UL < row_ + m_ ) )
      throw std::runtime_error( "Invalid transpose of an upper matrix" );

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   const ResultType tmp( trans(*this) );
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales all elements of the submatrix by the given scalar value \a scalar. Note
// that the function cannot be used to scale a submatrix on a lower or upper unitriangular matrix.
// The attempt to scale such a submatrix results in a compile time error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseSubmatrix<MT,aligned,true>& DenseSubmatrix<MT,aligned,true>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jend( column_ + n_ );

   for( size_t j=column_; j<jend; ++j )
   {
      const size_t ibegin( ( IsLower<MT>::value )
                           ?( ( IsStrictlyLower<MT>::value )
                              ?( max( j+1UL, row_ ) )
                              :( max( j, row_ ) ) )
                           :( row_ ) );
      const size_t iend  ( ( IsUpper<MT>::value )
                           ?( ( IsStrictlyUpper<MT>::value )
                              ?( min( j, row_+m_ ) )
                              :( min( j+1UL, row_+m_ ) ) )
                           :( row_+m_ ) );

      for( size_t i=ibegin; i<iend; ++i )
         matrix_(i,j) *= scalar;
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
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,aligned,true>::hasOverlap() const
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric<MT>::value, "Unsymmetric matrix detected" );

   if( ( row_ + m_ <= column_ ) || ( column_ + n_ <= row_ ) )
      return false;
   else return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying matrix.
//
// \param lhs The matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying matrix of type \a MT would be
// violated by an assignment of the given matrix \a rhs. In case the matrix would be preserved,
// the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the left-hand side dense matrix
        , bool SO2       // Storage order of the left-hand side dense matrix
        , typename MT3   // Type of the right-hand side matrix
        , bool SO3 >     // Storage order of the right-hand side matrix
inline typename EnableIf< Not< IsAdaptor<MT2> >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   UNUSED_PARAMETER( lhs, rhs );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying symmetric matrix.
//
// \param lhs The symmetric matrix to be assigned to.
// \param rhs The matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying symmetric matrix of type \a MT would
// be violated by an assignment of the given row-major matrix \a rhs. In case the matrix would be
// preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Type of the left-hand side dense matrix
        , bool SO2       // Storage order of the left-hand side dense matrix
        , typename MT3   // Type of the right-hand side matrix
        , bool SO3 >     // Storage order of the right-hand side matrix
inline typename EnableIf< And< IsSymmetric<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const Matrix<MT3,SO3>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( !hasOverlap() )
      return true;

   const bool   lower( row_ > column_ );
   const size_t size ( min( row_ + m_, column_ + n_ ) - ( lower ? row_ : column_ ) );

   if( size < 2UL )
      return true;

   const size_t row   ( lower ? 0UL : column_ - row_ );
   const size_t column( lower ? row_ - column_ : 0UL );

   return isSymmetric( submatrix( ~rhs, row, column, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;

      const size_t jbegin( ( containsDiagonal )
                           ?( ( IsStrictlyLower<MT2>::value )
                              ?( row_ + i - column_ )
                              :( row_ + i - column_ + 1UL ) )
                           :( 0UL ) );

      for( size_t j=jbegin; j<n_; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t ipos( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                         ?( column_ + j - row_ + 1UL )
                         :( column_ + j - row_ ) );
      const size_t iend( min( ipos, m_ ) );

      for( size_t i=0UL; i<iend; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniLower<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t ipos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( column_ + n_ - row_ )
                      :( column_ + n_ - row_ - 1UL ) );
   const size_t iend( min( ipos, m_ ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row_+i >= column_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                             ?( row_ + i - column_ )
                             :( row_ + i - column_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(i) );
      RhsIterator element( (~rhs).lowerBound( i, index ) );

      if( IsUniLower<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != row_+i-column_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying lower triangular matrix.
//
// \param lhs The lower matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying lower triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsLower<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( row_ + 1UL >= column_ + n_ )
      return true;

   const size_t jpos( ( IsUniLower<MT2>::value || IsStrictlyLower<MT2>::value )
                      ?( row_ - column_ )
                      :( row_ - column_ + 1UL ) );
   const size_t jbegin( ( row_ < column_ )?( 0UL ):( jpos ) );

   for( size_t j=jbegin; j<n_; ++j )
   {
      const bool containsDiagonal( column_+j < row_+m_ );

      const size_t index( ( IsStrictlyLower<MT2>::value && containsDiagonal )
                          ?( column_ + j - row_ + 1UL )
                          :( column_ + j - row_ ) );
      const RhsIterator last( (~rhs).lowerBound( min( index, m_ ), j ) );

      if( IsUniLower<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(j) || ( last->index() != column_+j-row_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(j); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t jpos( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                         ?( row_ + i - column_ + 1UL )
                         :( row_ + i - column_ ) );
      const size_t jend( min( jpos, n_ ) );

      for( size_t j=0UL; j<jend; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(i,row_+i-column_) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major dense matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      if( IsUniUpper<MT2>::value && containsDiagonal && !isOne( (~rhs)(column_+j-row_,j) ) )
         return false;

      const size_t ibegin( ( containsDiagonal )
                           ?( ( IsStrictlyUpper<MT2>::value )
                              ?( column_ + j - row_ )
                              :( column_ + j - row_ + 1UL ) )
                           :( 0UL ) );

      for( size_t i=ibegin; i<m_; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given row-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t ipos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( column_ - row_ )
                      :( column_ - row_ + 1UL ) );
   const size_t ibegin( ( column_ < row_ )?( 0UL ):( ipos ) );

   for( size_t i=ibegin; i<m_; ++i )
   {
      const bool containsDiagonal( row_+i < column_+n_ );

      const size_t index( ( IsStrictlyUpper<MT2>::value && containsDiagonal )
                          ?( row_ + i - column_ + 1UL )
                          :( row_ + i - column_ ) );
      const RhsIterator last( (~rhs).lowerBound( i, min( index, n_ ) ) );

      if( IsUniUpper<MT2>::value && containsDiagonal &&
          ( last == (~rhs).end(i) || ( last->index() != row_+i-column_ ) || !isOne( last->value() ) ) ) {
         return false;
      }

      for( RhsIterator element=(~rhs).begin(i); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying upper triangular matrix.
//
// \param lhs The upper matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying upper triangular matrix of type
// \a MT would be violated by an assignment of the given column-major sparse matrix \a rhs. In
// case the matrix would be preserved, the function returns \a true. Otherwise it returns
// \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< And< IsUpper<MT2>, Not< IsDiagonal<MT2> > >, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   if( column_ + 1UL >= row_ + m_ )
      return true;

   const size_t jpos( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                      ?( row_ + m_ - column_ )
                      :( row_ + m_ - column_ - 1UL ) );
   const size_t jend( min( jpos, n_ ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column_+j >= row_ );

      const size_t index( ( containsDiagonal )
                          ?( ( IsUniUpper<MT2>::value || IsStrictlyUpper<MT2>::value )
                             ?( column_ + j - row_ )
                             :( column_ + j - row_ + 1UL ) )
                          :( 0UL ) );

      const RhsIterator last( (~rhs).end(j) );
      RhsIterator element( (~rhs).lowerBound( index, j ) );

      if( IsUniUpper<MT2>::value && containsDiagonal ) {
         if( element == last || ( element->index() != column_+j-row_ ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<n_; ++j ) {
         if( ( row_ + i != column_ + j ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The dense matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major dense matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side dense matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const DenseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<m_; ++i ) {
         if( ( column_ + j != row_ + i ) && !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given row-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t i=0UL; i<m_; ++i ) {
      for( RhsIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element ) {
         if( ( row_ + i != column_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking for possible invariant violations of the underlying diagonal matrix.
//
// \param lhs The diagonal matrix to be assigned to.
// \param rhs The sparse matrix to be checked.
// \return \a true in case the invariants of the matrix are preserved, \a false if not.
//
// This function checks if the invariants of the underlying diagonal matrix of type \a MT would
// be violated by an assignment of the given column-major sparse matrix \a rhs. In case the matrix
// would be preserved, the function returns \a true. Otherwise it returns \a false.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2    // Type of the left-hand side dense matrix
        , bool SO2        // Storage order of the left-hand side dense matrix
        , typename MT3 >  // Type of the right-hand side sparse matrix
inline typename EnableIf< IsDiagonal<MT2>, bool >::Type
   DenseSubmatrix<MT,aligned,true>::preservesInvariant( const DenseMatrix<MT2,SO2>& lhs, const SparseMatrix<MT3,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT3 );

   UNUSED_PARAMETER( lhs );

   typedef typename MT3::ConstIterator  RhsIterator;

   for( size_t j=0UL; j<n_; ++j ) {
      for( RhsIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element ) {
         if( ( column_ + j != row_ + element->index() ) && !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
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
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,aligned,true>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense submatrix
        , bool AF2       // Alignment flag of the foreign dense submatrix
        , bool SO2 >     // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,aligned,true>::canAlias( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
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
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubmatrix<MT,aligned,true>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense submatrix
        , bool AF2       // Alignment flag of the foreign dense submatrix
        , bool SO2 >     // Storage order of the foreign dense submatrix
inline bool DenseSubmatrix<MT,aligned,true>::isAliased( const DenseSubmatrix<MT2,AF2,SO2>* alias ) const
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row_    + m_ > alias->row_    ) && ( row_    < alias->row_    + alias->m_ ) &&
            ( column_ + n_ > alias->column_ ) && ( column_ < alias->column_ + alias->n_ ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each column of the submatrix are guaranteed to
// conform to the alignment restrictions of the underlying element type.
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,aligned,true>::isAligned() const
{
   return true;
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
// rows and/or columns of the submatrix).
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseSubmatrix<MT,aligned,true>::canSMPAssign() const
{
   return ( columns() > SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the intrinsic element. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseSubmatrix<MT,aligned,true>::IntrinsicType
   DenseSubmatrix<MT,aligned,true>::load( size_t i, size_t j ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   return matrix_.load( row_+i, column_+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the intrinsic element. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseSubmatrix<MT,aligned,true>::IntrinsicType
   DenseSubmatrix<MT,aligned,true>::loadu( size_t i, size_t j ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   return matrix_.loadu( row_+i, column_+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store of a specific intrinsic element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of values
// inside the intrinsic element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void
   DenseSubmatrix<MT,aligned,true>::store( size_t i, size_t j, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   matrix_.store( row_+i, column_+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store of a specific intrinsic element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the row index must be a multiple of
// the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void
   DenseSubmatrix<MT,aligned,true>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   matrix_.storeu( row_+i, column_+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of an intrinsic element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific intrinsic element of
// the dense submatrix. The row index must be smaller than the number of rows and the column
// index must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the intrinsic element. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void
   DenseSubmatrix<MT,aligned,true>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows()         , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j < columns()      , "Invalid column access index" );

   matrix_.stream( row_+i, column_+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,aligned,true>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,true>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row_+i    ,column_+j) = (~rhs)(i    ,j);
         matrix_(row_+i+1UL,column_+j) = (~rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         matrix_(row_+ipos,column_+j) = (~rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,aligned,true>::BLAZE_TEMPLATE VectorizedAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,true>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   if( useStreaming &&
       m_*n_ > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(~rhs).isAliased( &matrix_ ) )
   {
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=0UL; i<m_; i+=IT::size )
            stream( i, j, (~rhs).load(i,j) );
   }
   else
   {
      const size_t ipos( m_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == ipos, "Invalid end calculation" );

      for( size_t j=0UL; j<n_; ++j ) {
         typename MT2::ConstIterator it( (~rhs).begin(j) );
         for( size_t i=0UL; i<ipos; i+=IT::size*4UL ) {
            store( i             , j, it.load() ); it += IT::size;
            store( i+IT::size    , j, it.load() ); it += IT::size;
            store( i+IT::size*2UL, j, it.load() ); it += IT::size;
            store( i+IT::size*3UL, j, it.load() ); it += IT::size;
         }
         for( size_t i=ipos; i<m_; i+=IT::size, it+=IT::size ) {
            store( i, j, it.load() );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,aligned,true>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row_+i,column_+j) = (~rhs)(i,j);
            }
         }
      }
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
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,true>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) = element->value();
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
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,true>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,aligned,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,true>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row_+i    ,column_+j) += (~rhs)(i    ,j);
         matrix_(row_+i+1UL,column_+j) += (~rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         matrix_(row_+ipos,column_+j) += (~rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,aligned,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,true>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      typename MT2::ConstIterator it( (~rhs).begin(j) );
      for( size_t i=0UL; i<ipos; i+=IT::size*4UL ) {
         store( i             , j, load(i             ,j) + it.load() ); it += IT::size;
         store( i+IT::size    , j, load(i+IT::size    ,j) + it.load() ); it += IT::size;
         store( i+IT::size*2UL, j, load(i+IT::size*2UL,j) + it.load() ); it += IT::size;
         store( i+IT::size*3UL, j, load(i+IT::size*3UL,j) + it.load() ); it += IT::size;
      }
      for( size_t i=ipos; i<m_; i+=IT::size, it+=IT::size ) {
         store( i, j, load(i,j) + it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,aligned,true>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row_+i,column_+j) += (~rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,true>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,true>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DenseSubmatrix<MT,aligned,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,true>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row_+i    ,column_+j) -= (~rhs)(i    ,j);
         matrix_(row_+i+1UL,column_+j) -= (~rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         matrix_(row_+ipos,column_+j) -= (~rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DenseSubmatrix<MT,aligned,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT2> >::Type
   DenseSubmatrix<MT,aligned,true>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t ipos( m_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == ipos, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      typename MT2::ConstIterator it( (~rhs).begin(j) );
      for( size_t i=0UL; i<ipos; i+=IT::size*4UL ) {
         store( i             , j, load(i             ,j) - it.load() ); it += IT::size;
         store( i+IT::size    , j, load(i+IT::size    ,j) - it.load() ); it += IT::size;
         store( i+IT::size*2UL, j, load(i+IT::size*2UL,j) - it.load() ); it += IT::size;
         store( i+IT::size*3UL, j, load(i+IT::size*3UL,j) - it.load() ); it += IT::size;
      }
      for( size_t i=ipos; i<m_; i+=IT::size, it+=IT::size ) {
         store( i, j, load(i,j) - it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side dense matrix
inline void DenseSubmatrix<MT,aligned,true>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row_+i,column_+j) -= (~rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,true>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         matrix_(row_+element->index(),column_+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >   // Type of the dense matrix
template< typename MT2 >  // Type of the right-hand side sparse matrix
inline void DenseSubmatrix<MT,aligned,true>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_RESTRICTED( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         matrix_(row_+i,column_+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  DENSESUBMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseSubmatrix operators */
//@{
template< typename MT, bool AF, bool SO >
inline void reset( DenseSubmatrix<MT,AF,SO>& dm );

template< typename MT, bool AF, bool SO >
inline void reset( DenseSubmatrix<MT,AF,SO>& dm, size_t i );

template< typename MT, bool AF, bool SO >
inline void clear( DenseSubmatrix<MT,AF,SO>& dm );

template< typename MT, bool AF, bool SO >
inline bool isDefault( const DenseSubmatrix<MT,AF,SO>& dm );

template< typename MT, bool AF, bool SO >
inline bool isSymmetric( const DenseSubmatrix<MT,AF,SO>& dm );

template< typename MT, bool AF, bool SO >
inline bool isLower( const DenseSubmatrix<MT,AF,SO>& dm );

template< typename MT, bool AF, bool SO >
inline bool isUpper( const DenseSubmatrix<MT,AF,SO>& dm );

template< typename MT, bool AF, bool SO >
inline bool isSame( const DenseSubmatrix<MT,AF,SO>& a, const DenseMatrix<MT,SO>& b );

template< typename MT, bool AF, bool SO >
inline bool isSame( const DenseMatrix<MT,SO>& a, const DenseSubmatrix<MT,AF,SO>& b );

template< typename MT, bool AF, bool SO >
inline bool isSame( const DenseSubmatrix<MT,AF,SO>& a, const DenseSubmatrix<MT,AF,SO>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dense submatrix.
// \ingroup dense_submatrix
//
// \param dm The dense submatrix to be resetted.
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void reset( DenseSubmatrix<MT,AF,SO>& dm )
{
   dm.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given dense submatrix.
// \ingroup dense_submatrix
//
// \param dm The dense submatrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given dense submatrix to
// their default value. In case the given submatrix is a \a rowMajor matrix the function resets
// the values in row \a i, if it is a \a columnMajor matrix the function resets the values in
// column \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void reset( DenseSubmatrix<MT,AF,SO>& dm, size_t i )
{
   dm.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dense matrix.
// \ingroup dense_submatrix
//
// \param dm The dense matrix to be cleared.
// \return void
//
// Clearing a dense submatrix is equivalent to resetting it via the reset() function.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline void clear( DenseSubmatrix<MT,AF,SO>& dm )
{
   dm.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense submatrix is in default state.
// \ingroup dense_submatrix
//
// \param dm The dense submatrix to be tested for its default state.
// \return \a true in case the given submatrix is component-wise zero, \a false otherwise.
//
// This function checks whether the submatrix is in default state. For instance, in
// case the submatrix is instantiated for a built-in integral or floating point data type, the
// function returns \a true in case all submatrix elements are 0 and \a false in case any submatrix
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isDefault( const DenseSubmatrix<MT,AF,SO>& dm )
{
   using blaze::isDefault;

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<(~dm).rows(); ++i )
         for( size_t j=0UL; j<(~dm).columns(); ++j )
            if( !isDefault( (~dm)(i,j) ) )
               return false;
   }
   else {
      for( size_t j=0UL; j<(~dm).columns(); ++j )
         for( size_t i=0UL; i<(~dm).rows(); ++i )
            if( !isDefault( (~dm)(i,j) ) )
               return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense submatrix is symmetric.
// \ingroup dense_submatrix
//
// \param dm The dense submatrix to be checked.
// \return \a true if the submatrix is symmetric, \a false if not.
//
// This function checks if the given dense submatrix is symmetric. The submatrix is considered to
// be symmetric if it is a square matrix whose transpose is equal to itself (\f$ A = A^T \f$). The
// following code example demonstrates the use of the function:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  Matrix;

   Matrix A( 32UL, 16UL );
   // ... Initialization

   blaze::DenseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isSymmetric( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSymmetric( const DenseSubmatrix<MT,AF,SO>& dm )
{
   typedef DenseMatrix< DenseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsSymmetric<MT>::value && dm.row_ == dm.column_ && dm.m_ == dm.n_ )
      return true;
   else return isSymmetric( static_cast<const BaseType&>( dm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense submatrix is a lower triangular matrix.
// \ingroup dense_submatrix
//
// \param dm The dense submatrix to be checked.
// \return \a true if the submatrix is a lower triangular matrix, \a false if not.
//
// This function checks if the given dense submatrix is a lower triangular matrix. The matrix is
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

   blaze::DenseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isLower( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isLower( const DenseSubmatrix<MT,AF,SO>& dm )
{
   typedef DenseMatrix< DenseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsLower<MT>::value && dm.row_ == dm.column_ && dm.m_ == dm.n_ )
      return true;
   else return isLower( static_cast<const BaseType&>( dm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense submatrix is an upper triangular matrix.
// \ingroup dense_submatrix
//
// \param dm The dense submatrix to be checked.
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

   blaze::DenseSubmatrix<Matrix> sm( A, 8UL, 8UL, 16UL, 16UL );

   if( isUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isUpper( const DenseSubmatrix<MT,AF,SO>& dm )
{
   typedef DenseMatrix< DenseSubmatrix<MT,AF,SO>, SO >  BaseType;

   if( IsUpper<MT>::value && dm.row_ == dm.column_ && dm.m_ == dm.n_ )
      return true;
   else return isUpper( static_cast<const BaseType&>( dm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense matrix and submatrix represent the same observable state.
// \ingroup dense_submatrix
//
// \param a The dense submatrix to be tested for its state.
// \param b The dense matrix to be tested for its state.
// \return \a true in case the dense submatrix and matrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given submatrix refers to the full given
// dense matrix and by that represents the same observable state. In this case, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSame( const DenseSubmatrix<MT,AF,SO>& a, const DenseMatrix<MT,SO>& b )
{
   return ( isSame( a.matrix_, ~b ) && ( a.rows() == (~b).rows() ) && ( a.columns() == (~b).columns() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense matrix and submatrix represent the same observable state.
// \ingroup dense_submatrix
//
// \param a The dense matrix to be tested for its state.
// \param b The dense submatrix to be tested for its state.
// \return \a true in case the dense matrix and submatrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given submatrix refers to the full given
// dense matrix and by that represents the same observable state. In this case, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSame( const DenseMatrix<MT,SO>& a, const DenseSubmatrix<MT,AF,SO>& b )
{
   return ( isSame( ~a, b.matrix_ ) && ( (~a).rows() == b.rows() ) && ( (~a).columns() == b.columns() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given submatrices represent the same observable state.
// \ingroup dense_submatrix
//
// \param a The first dense submatrix to be tested for its state.
// \param b The second dense submatrix to be tested for its state.
// \return \a true in case the two submatrices share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given submatrices refer to exactly the
// same part of the same dense matrix. In case both submatrices represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline bool isSame( const DenseSubmatrix<MT,AF,SO>& a, const DenseSubmatrix<MT,AF,SO>& b )
{
   return ( isSame( a.matrix_, b.matrix_ ) &&
            ( a.row_ == b.row_ ) && ( a.column_ == b.column_ ) &&
            ( a.m_ == b.m_ ) && ( a.n_ == b.n_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given dense submatrix.
// \ingroup dense_submatrix
//
// \param dm The submatrix to be derestricted.
// \return Submatrix without access restrictions.
//
// This function removes all restrictions on the data access to the given submatrix. It returns a
// submatrix that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool AF      // Alignment flag
        , bool SO >    // Storage order
inline typename DerestrictTrait< DenseSubmatrix<MT,AF,SO> >::Type
   derestrict( DenseSubmatrix<MT,AF,SO>& dm )
{
   typedef typename DerestrictTrait< DenseSubmatrix<MT,AF,SO> >::Type  ReturnType;
   return ReturnType( derestrict( dm.matrix_ ), dm.row_, dm.column_, dm.m_, dm.n_ );
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
/*!\brief Creating a view on a specific submatrix of another dense submatrix.
// \ingroup views
//
// \param dm The constant dense submatrix
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the other dense submatrix.
//
// This function returns an expression representing the specified submatrix of the given
// dense submatrix.
*/
template< bool AF1     // Required alignment flag
        , typename MT  // Type of the sparse submatrix
        , bool AF2     // Present alignment flag
        , bool SO >    // Storage order
inline const DenseSubmatrix<MT,AF1,SO>
   submatrix( const DenseSubmatrix<MT,AF2,SO>& dm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   if( ( row + m > dm.rows() ) || ( column + n > dm.columns() ) )
      throw std::invalid_argument( "Invalid submatrix specification" );

   return DenseSubmatrix<MT,AF1,SO>( dm.matrix_, dm.row_ + row, dm.column_ + column, m, n );
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
struct IsRestricted< DenseSubmatrix<MT,AF,SO> > : public If< IsRestricted<MT>, TrueType, FalseType >::Type
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
struct DerestrictTrait< DenseSubmatrix<MT,AF,SO> >
{
   typedef DenseSubmatrix< typename RemoveReference< typename DerestrictTrait<MT>::Type >::Type, AF >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF, bool SO >
struct HasConstDataAccess< DenseSubmatrix<MT,AF,SO> >
   : public If< HasConstDataAccess<MT>, TrueType, FalseType >::Type
{
   enum { value = HasConstDataAccess<MT>::value };
   typedef typename If< HasConstDataAccess<MT>, TrueType, FalseType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF, bool SO >
struct HasMutableDataAccess< DenseSubmatrix<MT,AF,SO> >
   : public If< HasMutableDataAccess<MT>, TrueType, FalseType >::Type
{
   enum { value = HasMutableDataAccess<MT>::value };
   typedef typename If< HasMutableDataAccess<MT>, TrueType, FalseType >::Type  Type;
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
struct SubmatrixTrait< DenseSubmatrix<MT,AF,SO> >
{
   typedef typename SubmatrixTrait< typename DenseSubmatrix<MT,AF,SO>::ResultType >::Type  Type;
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
struct SubmatrixExprTrait< DenseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef DenseSubmatrix<MT,AF2,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< const DenseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef DenseSubmatrix<MT,AF2,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< volatile DenseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef DenseSubmatrix<MT,AF2,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool AF1, bool SO, bool AF2 >
struct SubmatrixExprTrait< const volatile DenseSubmatrix<MT,AF1,SO>, AF2 >
{
   typedef DenseSubmatrix<MT,AF2,SO>  Type;
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
struct RowTrait< DenseSubmatrix<MT,AF,SO> >
{
   typedef typename RowTrait< typename DenseSubmatrix<MT,AF,SO>::ResultType >::Type  Type;
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
struct ColumnTrait< DenseSubmatrix<MT,AF,SO> >
{
   typedef typename ColumnTrait< typename DenseSubmatrix<MT,AF,SO>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
