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
#include <stdexcept>
#include <vector>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/Submatrix.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/AlignmentFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>


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
   typedef typename SelectType< IsExpression<MT>::value, MT, MT& >::Type  Operand;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given sparse matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the sparse submatrix will
       return references and iterators to const. Otherwise \a useConst will be set to 0 and
       the sparse submatrix will offer write access to the sparse matrix elements both via
       the function call operator and iterators. */
   enum { useConst = IsConst<MT>::value };
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
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;
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

    public:
      //**Type definitions*************************************************************************
      //! Return type of the value member function.
      typedef typename SelectType< returnConst, const ElementType&, ElementType& >::Type  ReferenceType;
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
      inline ReferenceType value() const {
         return pos_->value();
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
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
   typedef typename SelectType< useConst, ConstIterator, SubmatrixIterator<MT,typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n );
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
                                      inline SparseSubmatrix& operator= ( const SparseSubmatrix& rhs );
   template< typename MT2, bool SO2 > inline SparseSubmatrix& operator= ( const DenseMatrix<MT2,SO2>&  rhs );
   template< typename MT2, bool SO2 > inline SparseSubmatrix& operator= ( const SparseMatrix<MT2,SO2>& rhs );
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
                              inline size_t           rows() const;
                              inline size_t           columns() const;
                              inline size_t           capacity() const;
                              inline size_t           capacity( size_t i ) const;
                              inline size_t           nonZeros() const;
                              inline size_t           nonZeros( size_t i ) const;
                              inline void             reset();
                              inline void             reset( size_t i );
                                     Iterator         insert( size_t i, size_t j, const ElementType& value );
                              inline void             erase( size_t i, size_t j );
                              inline Iterator         erase( size_t i, Iterator pos );
                              inline Iterator         erase( size_t i, Iterator first, Iterator last );
                              inline void             reserve( size_t nonzeros );
                                     void             reserve( size_t i, size_t nonzeros );
                              inline void             trim();
                              inline void             trim( size_t i );
   template< typename Other > inline SparseSubmatrix& scale( Other scalar );
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
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
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
// \param row The index of the first row of the submatrix in the given sparse matrix.
// \param column The index of the first column of the submatrix in the given sparse matrix.
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
inline SparseSubmatrix<MT,AF,SO>::SparseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_( matrix )  // The sparse matrix containing the submatrix
   , row_   ( row    )  // The first row of the submatrix
   , column_( column )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
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
/*!\brief 2D-access to the sparse submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
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
//
// The sparse submatrix is initialized as a copy of the given sparse submatrix. In case the
// current sizes of the two submatrices don't match, a \a std::invalid_argument exception is
// thrown.
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

   if( rows() != rhs.rows() || columns() != rhs.columns() )
      throw std::invalid_argument( "Submatrix sizes do not match" );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
      assign( *this, rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for dense matrices.
//
// \param rhs Dense matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The sparse submatrix is initialized as a copy of the given dense matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator=( const DenseMatrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( typename MT2::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( RequiresEvaluation<MT2>::value || (~rhs).canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( ~rhs );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different sparse matrices.
//
// \param rhs Sparse matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The sparse submatrix is initialized as a copy of the given sparse matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand side sparse matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator=( const SparseMatrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( typename MT2::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( RequiresEvaluation<MT2>::value || (~rhs).canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( ~rhs );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator+=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::addAssign;

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   addAssign( *this, ~rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator-=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::subAssign;

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   subAssign( *this, ~rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the sparse submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF       // Alignment flag
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,SO>&
   SparseSubmatrix<MT,AF,SO>::operator*=( const Matrix<MT2,SO2>& rhs )
{
   if( columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE        ( MultType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType   );

   const MultType tmp( *this * (~rhs) );
   reset();
   assign( tmp );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a sparse submatrix
//        and a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the sparse submatrix.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,SO> >::Type&
   SparseSubmatrix<MT,AF,SO>::operator*=( Other rhs )
{
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
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,SO> >::Type&
   SparseSubmatrix<MT,AF,SO>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
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
   for( size_t i=row_; i<row_+m_; ++i ) {
      matrix_.erase( i, matrix_.lowerBound( i, column_ ), matrix_.lowerBound( i, column_+n_ ) );
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
   matrix_.erase( index, matrix_.lowerBound( index, column_ ), matrix_.lowerBound( index, column_+n_ ) );
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
typename SparseSubmatrix<MT,AF,SO>::Iterator
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
/*!\brief Scaling of the sparse submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the sparse submatrix.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline SparseSubmatrix<MT,AF,SO>& SparseSubmatrix<MT,AF,SO>::scale( Other scalar )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
         element->value() *= scalar;
   }

   return *this;
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

   B.reserve( 3 );       // Reserving enough space for 3 non-zero elements
   B.append( 0, 1, 1 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );      // Finalizing row 0
   B.append( 1, 1, 2 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );      // Finalizing row 1
   B.append( 2, 0, 3 );  // Appending the value 3 in row 2 with column index 0
   B.finalize( 2 );      // Finalizing row 2
   \endcode

// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
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
// \b Note: Although finalize() does not allocate new memory, it still invalidates all iterators
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
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the conAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF           // Alignment flag
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,AF,SO>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
         append( element->index(), j, element->value() );
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

   const AddType tmp( *this + (~rhs) );
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

   const AddType tmp( *this + (~rhs) );
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

   const SubType tmp( *this - (~rhs) );
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

   const SubType tmp( *this - (~rhs) );
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
   typedef typename SelectType< IsExpression<MT>::value, MT, MT& >::Type  Operand;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given sparse matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the sparse submatrix will
       return references and iterators to const. Otherwise \a useConst will be set to 0 and
       the sparse submatrix will offer write access to the sparse matrix elements both via
       the function call operator and iterators. */
   enum { useConst = IsConst<MT>::value };
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
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;
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

    public:
      //**Type definitions*************************************************************************
      typedef typename SelectType< returnConst, const ElementType&, ElementType& >::Type  ReferenceType;
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
      inline ReferenceType value() const {
         return pos_->value();
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
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
   typedef typename SelectType< useConst, ConstIterator, SubmatrixIterator<MT,typename MT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n );
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
                                     inline SparseSubmatrix& operator= ( const SparseSubmatrix& rhs );
   template< typename MT2, bool SO > inline SparseSubmatrix& operator= ( const DenseMatrix<MT2,SO>&  rhs );
   template< typename MT2, bool SO > inline SparseSubmatrix& operator= ( const SparseMatrix<MT2,SO>& rhs );
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
                              inline size_t           rows() const;
                              inline size_t           columns() const;
                              inline size_t           capacity() const;
                              inline size_t           capacity( size_t i ) const;
                              inline size_t           nonZeros() const;
                              inline size_t           nonZeros( size_t i ) const;
                              inline void             reset();
                              inline void             reset( size_t i );
                                     Iterator         insert( size_t i, size_t j, const ElementType& value );
                              inline void             erase( size_t i, size_t j );
                              inline Iterator         erase( size_t i, Iterator pos );
                              inline Iterator         erase( size_t i, Iterator first, Iterator last );
                              inline void             reserve( size_t nonzeros );
                                     void             reserve( size_t i, size_t nonzeros );
                              inline void             trim();
                              inline void             trim( size_t j );
   template< typename Other > inline SparseSubmatrix& scale( Other scalar );
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
   template< typename MT2, bool SO > inline void assign   ( const DenseMatrix<MT2,SO>&     rhs );
   template< typename MT2 >          inline void assign   ( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 >          inline void assign   ( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2, bool SO > inline void addAssign( const DenseMatrix<MT2,SO>&     rhs );
   template< typename MT2, bool SO > inline void addAssign( const SparseMatrix<MT2,SO>&    rhs );
   template< typename MT2, bool SO > inline void subAssign( const DenseMatrix<MT2,SO>&     rhs );
   template< typename MT2, bool SO > inline void subAssign( const SparseMatrix<MT2,SO>&    rhs );
   //@}
   //**********************************************************************************************

 private:
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
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
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
// \param row The index of the first row of the submatrix in the given sparse matrix.
// \param column The index of the first column of the submatrix in the given sparse matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// In case the submatrix is not properly specified (i.e. if the specified submatrix is not
// contained in the given sparse matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool AF >    // Alignment flag
inline SparseSubmatrix<MT,AF,true>::SparseSubmatrix( Operand matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_( matrix )  // The sparse matrix containing the submatrix
   , row_   ( row    )  // The first row of the submatrix
   , column_( column )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
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
/*!\brief 2D-access to the sparse submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
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
//
// The sparse submatrix is initialized as a copy of the given sparse submatrix. In case the
// current sizes of the two submatrices don't match, a \a std::invalid_argument exception is
// thrown.
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

   if( rows() != rhs.rows() || columns() != rhs.columns() )
      throw std::invalid_argument( "Submatrix sizes do not match" );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
      assign( *this, rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for dense matrices.
//
// \param rhs Dense matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The sparse submatrix is initialized as a copy of the given dense matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator=( const DenseMatrix<MT2,SO>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( RequiresEvaluation<MT2>::value || (~rhs).canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( ~rhs );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
      assign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different sparse matrices.
//
// \param rhs Sparse matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The sparse submatrix is initialized as a copy of the given sparse matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO >     // Storage order of the right-hand side sparse matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator=( const SparseMatrix<MT2,SO>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( typename MT2::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( RequiresEvaluation<MT2>::value || (~rhs).canAlias( &matrix_ ) ) {
      const typename MT2::ResultType tmp( ~rhs );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
      assign( *this, ~rhs );
   }

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
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator+=( const Matrix<MT2,SO>& rhs )
{
   using blaze::addAssign;

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   addAssign( *this, ~rhs );

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
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator-=( const Matrix<MT2,SO>& rhs )
{
   using blaze::subAssign;

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   subAssign( *this, ~rhs );

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
//
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool AF >     // Alignment flag
template< typename MT2  // Type of the right-hand side matrix
        , bool SO >     // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,AF,true>&
   SparseSubmatrix<MT,AF,true>::operator*=( const Matrix<MT2,SO>& rhs )
{
   if( columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   typedef typename MultTrait<ResultType,typename MT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE        ( MultType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType   );

   const MultType tmp( *this * (~rhs) );
   reset();
   assign( tmp );

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
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,true> >::Type&
   SparseSubmatrix<MT,AF,true>::operator*=( Other rhs )
{
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
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,AF,true> >::Type&
   SparseSubmatrix<MT,AF,true>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
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
   for( size_t j=column_; j<column_+n_; ++j ) {
      matrix_.erase( j, matrix_.lowerBound( row_, j ), matrix_.lowerBound( row_+m_, j ) );
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
   matrix_.erase( index, matrix_.lowerBound( row_, index ), matrix_.lowerBound( row_+m_, index ) );
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
typename SparseSubmatrix<MT,AF,true>::Iterator
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
/*!\brief Scaling of the sparse submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the sparse submatrix.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the scalar value
inline SparseSubmatrix<MT,AF,true>& SparseSubmatrix<MT,AF,true>::scale( Other scalar )
{
   for( size_t i=0UL; i<columns(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
         element->value() *= scalar;
   }

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

   B.reserve( 3 );       // Reserving enough space for 3 non-zero elements
   B.append( 0, 1, 1 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );      // Finalizing row 0
   B.append( 1, 1, 2 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );      // Finalizing row 1
   B.append( 2, 0, 3 );  // Appending the value 3 in row 2 with column index 0
   B.finalize( 2 );      // Finalizing row 2
   \endcode

// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
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
// \b Note: Although finalize() does not allocate new memory, it still invalidates all iterators
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
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
// to the conAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool AF >         // Alignment flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,AF,true>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
         append( i, j, (~rhs)(i,j), true );
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
         append( i, element->index(), element->value() );
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
         append( element->index(), j, element->value(), true );
      }
      finalize( j );
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

   const AddType tmp( *this + (~rhs) );
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

   const AddType tmp( *this + (~rhs) );
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

   const SubType tmp( *this - (~rhs) );
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

   const SubType tmp( *this - (~rhs) );
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
inline void clear( SparseSubmatrix<MT,AF,SO>& sm );

template< typename MT, bool AF, bool SO >
inline bool isDefault( const SparseSubmatrix<MT,AF,SO>& sm );
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

   if( ( row + m > sm.rows() ) || ( column + n > sm.columns() ) )
      throw std::invalid_argument( "Invalid submatrix specification" );

   return SparseSubmatrix<MT,AF1,SO>( sm.matrix_, sm.row_ + row, sm.column_ + column, m, n );
}
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
