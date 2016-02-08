//=================================================================================================
/*!
//  \file blaze/math/views/DenseRow.h
//  \brief Header file for the DenseRow class template
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

#ifndef _BLAZE_MATH_VIEWS_DENSEROW_H_
#define _BLAZE_MATH_VIEWS_DENSEROW_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Row.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Null.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_row DenseRow
// \ingroup views
*/
/*!\brief Reference to a specific row of a dense matrix.
// \ingroup dense_row
//
// The DenseRow template represents a reference to a specific row of a dense matrix primitive.
// The type of the dense matrix is specified via the first template parameter:

   \code
   template< typename MT, bool SO, bool SF >
   class DenseRow;
   \endcode

//  - MT: specifies the type of the dense matrix primitive. DenseRow can be used with every dense
//        matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the dense matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//  - SF: specifies whether the given matrix is a symmetric matrix or not. Also this parameter
//        doesn't have to be explicitly defined, but is automatically derived from the first
//        template parameter.
//
//
// \n \section dense_row_setup Setup of Dense Rows
//
// A reference to a dense row can be created very conveniently via the \c row() function. This
// reference can be treated as any other row vector, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The reference can also be used on
// both sides of an assignment: The row can either be used as an alias to grant write access to a
// specific row of a matrix primitive on the left-hand side of an assignment or to grant read-access
// to a specific row of a matrix primitive or expression on the right-hand side of an assignment.
// The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>      DenseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   DenseMatrixType  A, B;
   // ... Resizing and initialization

   // Setting the 2nd row of matrix A to x
   blaze::DenseRow<DenseMatrixType> row2 = row( A, 2UL );
   row2 = x;

   // Setting the 3rd row of matrix B to y
   row( B, 3UL ) = y;

   // Setting x to the 1st row of matrix B
   x = row( B, 1UL );

   // Setting y to the 4th row of the result of the matrix multiplication
   y = row( A * B, 4UL );
   \endcode

// \n \section dense_row_element_access Element access
//
// A dense row can be used like any other row vector. For instance, the elements of the dense
// row can be directly accessed with the subscript operator.

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the 4th row of matrix A
   blaze::DenseRow<MatrixType> row4 = row( A, 4UL );

   // Setting the 1st element of the dense row, which corresponds
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
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseRow<MatrixType>                RowType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of matrix A
   RowType row31 = row( A, 31UL );

   for( RowType::Iterator it=row31.begin(); it!=row31.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense row value
      ... = *it;  // OK: Read access to the dense row value.
   }

   for( RowType::ConstIterator it=row31.begin(); it!=row31.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense row value.
   }
   \endcode

// \n \section dense_row_common_operations Common Operations
//
// The current number of row elements can be obtained via the \c size() function, the current
// capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since rows are references to specific rows of a matrix,
// several operations are not possible on views, such as resizing and swapping:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseRow<MatrixType>                RowType;

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

// \n \section dense_row_arithmetic_operations Arithmetic Operations
//
// The following example gives an impression of the use of DenseRow within arithmetic operations.
// All operations (addition, subtraction, multiplication, scaling, ...) can be performed on all
// possible combinations of dense and sparse vectors with fitting element types:

   \code
   blaze::DynamicVector<double,blaze::rowVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::rowVector> c( 2UL );
   c[1] = 3.0;

   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrix;
   DenseMatrix A( 4UL, 2UL );  // Non-initialized 4x2 matrix

   typedef blaze::DenseRow<DenseMatrix>  RowType;
   RowType row0( row( A, 0UL ) );  // Reference to the 0th row of A

   row0[0] = 0.0;        // Manual initialization of the 0th row of A
   row0[1] = 0.0;
   row( A, 1UL ) = 1.0;  // Homogeneous initialization of the 1st row of A
   row( A, 2UL ) = a;    // Dense vector initialization of the 2nd row of A
   row( A, 3UL ) = c;    // Sparse vector initialization of the 3rd row of A

   b = row0 + a;              // Dense vector/dense vector addition
   b = c + row( A, 1UL );     // Sparse vector/dense vector addition
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

// \n \section dense_row_on_column_major_matrix Dense Row on a Column-Major Matrix
//
// It is especially noteworthy that row views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// also possible to traverse a row of a column-major matrix. For instance:

   \code
   typedef blaze::DynamicMatrix<int,blaze::columnMajor>  MatrixType;
   typedef blaze::DenseRow<MatrixType>                   RowType;

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
   blaze::DynamicMatrix<double,blaze::columnMajor> A( 128UL, 128UL );
   blaze::DynamicMatrix<double,blaze::columnMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th row of the multiplication between A and B ...
   blaze::DynamicVector<double,blaze::rowVector> x = row( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // the 15th row of the column-major matrix A with B.
   blaze::DynamicVector<double,blaze::rowVector> x = row( A, 15UL ) * B;
   \endcode

// Although Blaze performs the resulting vector/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix A would result in a more efficient evaluation.
*/
template< typename MT                            // Type of the dense matrix
        , bool SO = IsRowMajorMatrix<MT>::value  // Storage order
        , bool SF = IsSymmetric<MT>::value >     // Symmetry flag
class DenseRow : public DenseVector< DenseRow<MT,SO,SF>, true >
               , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the row element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseRow<MT,SO,SF>                  This;           //!< Type of this DenseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseRow&                     CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant row value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant row value.
   typedef typename If< Or< IsConst<MT>, Not< HasMutableDataAccess<MT> > >
                      , ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, typename MT::Iterator >::Type  Iterator;
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
   explicit inline DenseRow( MT& matrix, size_t index );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseRow& operator=( const ElementType& rhs );
   inline DenseRow& operator=( const DenseRow& rhs );

   template< typename VT > inline DenseRow& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator*=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline DenseRow& operator*=( const SparseVector<VT,true>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseRow >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseRow >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t    size() const;
                              inline size_t    capacity() const;
                              inline size_t    nonZeros() const;
                              inline void      reset();
   template< typename Other > inline DenseRow& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedAddAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedSubAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedMultAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::multiplication };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool SO2, bool SF2 >
   inline bool canAlias( const DenseRow<MT2,SO2,SF2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool SO2, bool SF2 >
   inline bool isAliased( const DenseRow<MT2,SO2,SF2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   BLAZE_ALWAYS_INLINE IntrinsicType load ( size_t index ) const;
   BLAZE_ALWAYS_INLINE IntrinsicType loada( size_t index ) const;
   BLAZE_ALWAYS_INLINE IntrinsicType loadu( size_t index ) const;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void storea( size_t index, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void stream( size_t index, const IntrinsicType& value );

   template< typename VT >
   inline typename DisableIf< VectorizedAssign<VT> >::Type
      assign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedAssign<VT> >::Type
      assign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void assign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedAddAssign<VT> >::Type
      addAssign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedAddAssign<VT> >::Type
      addAssign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void addAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedSubAssign<VT> >::Type
      subAssign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedSubAssign<VT> >::Type
      subAssign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void subAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedMultAssign<VT> >::Type
      multAssign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedMultAssign<VT> >::Type
      multAssign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void multAssign( const SparseVector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The dense matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename MT2, bool SO2, bool SF2 > friend class DenseRow;

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const DenseRow<MT2,SO2,SF2>& row );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const DenseRow<MT2,SO2,SF2>& a, const DenseRow<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< DenseRow<MT2,SO2,SF2> >::Type
      derestrict( DenseRow<MT2,SO2,SF2>& dm );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
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
/*!\brief The constructor for DenseRow.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline DenseRow<MT,SO,SF>::DenseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the row
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::Reference DenseRow<MT,SO,SF>::operator[]( size_t index )
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstReference
   DenseRow<MT,SO,SF>::operator[]( size_t index ) const
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::Reference DenseRow<MT,SO,SF>::at( size_t index )
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstReference DenseRow<MT,SO,SF>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::Pointer DenseRow<MT,SO,SF>::data()
{
   return matrix_.data( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstPointer DenseRow<MT,SO,SF>::data() const
{
   return matrix_.data( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::Iterator DenseRow<MT,SO,SF>::begin()
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstIterator DenseRow<MT,SO,SF>::begin() const
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstIterator DenseRow<MT,SO,SF>::cbegin() const
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::Iterator DenseRow<MT,SO,SF>::end()
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstIterator DenseRow<MT,SO,SF>::end() const
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DenseRow<MT,SO,SF>::ConstIterator DenseRow<MT,SO,SF>::cend() const
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
/*!\brief Homogenous assignment to all row elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row.
//
// This function homogeneously assigns the given value to all elements of the row. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator=( const ElementType& rhs )
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

   for( size_t j=jbegin; j<jend; ++j )
      matrix_(row_,j) = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DenseRow.
//
// \param rhs Dense row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator=( const DenseRow& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, rhs );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator=( const Vector<VT,true>& rhs )
{
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
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAddAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !trySubAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a dense vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator*=( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryMultAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::operator*=( const SparseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const ResultType right( *this * (~rhs) );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, right );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a dense row and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the vector.
//
// This operator cannot be used for rows on lower or upper unitriangular matrices. The attempt
// to scale such a row results in a compilation error!
*/
template< typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,SO,SF> >::Type&
   DenseRow<MT,SO,SF>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   return operator=( (*this) * rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense row by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the vector.
//
// This operator cannot be used for rows on lower or upper unitriangular matrices. The attempt
// to scale such a row results in a compilation error!
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,SO,SF> >::Type&
   DenseRow<MT,SO,SF>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   return operator=( (*this) / rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the row.
//
// \return The size of the row.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t DenseRow<MT,SO,SF>::size() const
{
   return matrix_.columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the dense row.
//
// \return The capacity of the dense row.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t DenseRow<MT,SO,SF>::capacity() const
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
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline size_t DenseRow<MT,SO,SF>::nonZeros() const
{
   return matrix_.nonZeros( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void DenseRow<MT,SO,SF>::reset()
{
   matrix_.reset( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the scalar value
inline DenseRow<MT,SO,SF>& DenseRow<MT,SO,SF>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jbegin( ( IsUpper<MT>::value )
                        ?( ( IsStrictlyUpper<MT>::value )
                           ?( row_+1UL )
                           :( row_ ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower<MT>::value )
                        ?( ( IsStrictlyLower<MT>::value )
                           ?( row_ )
                           :( row_+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_(row_,j) *= scalar;
   }

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the dense row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,SO,SF>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense row can alias with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT   // Type of the dense matrix
        , bool SO       // Storage order
        , bool SF >     // Symmetry flag
template< typename MT2  // Data type of the foreign dense row
        , bool SO2      // Storage order of the foreign dense row
        , bool SF2 >    // Symmetry flag of the foreign dense row
inline bool DenseRow<MT,SO,SF>::canAlias( const DenseRow<MT2,SO2,SF2>* alias ) const
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row_ == alias->row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF >         // Symmetry flag
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,SO,SF>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense row is aliased with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT   // Type of the dense matrix
        , bool SO       // Storage order
        , bool SF >     // Symmetry flag
template< typename MT2  // Data type of the foreign dense row
        , bool SO2      // Storage order of the foreign dense row
        , bool SF2 >    // Symmetry flag of the foreign dense row
inline bool DenseRow<MT,SO,SF>::isAliased( const DenseRow<MT2,SO2,SF2>* alias ) const
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row_ == alias->row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense row is properly aligned in memory.
//
// \return \a true in case the dense row is aligned, \a false if not.
//
// This function returns whether the dense row is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense row are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool DenseRow<MT,SO,SF>::isAligned() const
{
   return matrix_.isAligned();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense row can be used in SMP assignments.
//
// \return \a true in case the dense row can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense row).
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool DenseRow<MT,SO,SF>::canSMPAssign() const
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Load of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded intrinsic element.
//
// This function performs a load of a specific intrinsic element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE typename DenseRow<MT,SO,SF>::IntrinsicType
   DenseRow<MT,SO,SF>::load( size_t index ) const
{
   return matrix_.load( row_, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE typename DenseRow<MT,SO,SF>::IntrinsicType
   DenseRow<MT,SO,SF>::loada( size_t index ) const
{
   return matrix_.loada( row_, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE typename DenseRow<MT,SO,SF>::IntrinsicType
   DenseRow<MT,SO,SF>::loadu( size_t index ) const
{
   return matrix_.loadu( row_, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs a store a specific intrinsic element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE void DenseRow<MT,SO,SF>::store( size_t index, const IntrinsicType& value )
{
   matrix_.store( row_, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE void DenseRow<MT,SO,SF>::storea( size_t index, const IntrinsicType& value )
{
   matrix_.storea( row_, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unligned store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE void DenseRow<MT,SO,SF>::storeu( size_t index, const IntrinsicType& value )
{
   matrix_.storeu( row_, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific intrinsic element of the
// dense row. The index must be smaller than the number of matrix columns. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized evaluation
// of expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
BLAZE_ALWAYS_INLINE void DenseRow<MT,SO,SF>::stream( size_t index, const IntrinsicType& value )
{
   matrix_.stream( row_, index, value );
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseRow<MT,SO,SF>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) = (~rhs)[j    ];
      matrix_(row_,j+1UL) = (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) = (~rhs)[jpos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseRow<MT,SO,SF>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t columns( size() );

   const size_t jpos( ( remainder )?( columns & size_t(-IT::size) ):( columns ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( columns - ( columns % (IT::size) ) ) == jpos, "Invalid end calculation" );

   if( useStreaming && columns > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(~rhs).isAliased( &matrix_ ) )
   {
      size_t j( 0UL );

      for( ; j<jpos; j+=IT::size ) {
         matrix_.stream( row_, j, (~rhs).load(j) );
      }
      for( ; remainder && j<columns; ++j ) {
         matrix_(row_,j) = (~rhs)[j];
      }
   }
   else
   {
      size_t j( 0UL );
      typename VT::ConstIterator it( (~rhs).begin() );

      for( ; (j+IT::size*3UL) < jpos; j+=IT::size*4UL ) {
         matrix_.store( row_, j             , it.load() ); it += IT::size;
         matrix_.store( row_, j+IT::size    , it.load() ); it += IT::size;
         matrix_.store( row_, j+IT::size*2UL, it.load() ); it += IT::size;
         matrix_.store( row_, j+IT::size*3UL, it.load() ); it += IT::size;
      }
      for( ; j<jpos; j+=IT::size, it+=IT::size ) {
         matrix_.store( row_, j, it.load() );
      }
      for( ; remainder && j<columns; ++j, ++it ) {
         matrix_(row_,j) = *it;
      }
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO,SF>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) = element->value();
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseRow<MT,SO,SF>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) += (~rhs)[j    ];
      matrix_(row_,j+1UL) += (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) += (~rhs)[jpos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseRow<MT,SO,SF>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t columns( size() );

   const size_t jpos( ( remainder )?( columns & size_t(-IT::size) ):( columns ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( columns - ( columns % (IT::size) ) ) == jpos, "Invalid end calculation" );

   size_t j( 0UL );
   typename VT::ConstIterator it( (~rhs).begin() );

   for( ; (j+IT::size*3UL) < jpos; j+=IT::size*4UL ) {
      matrix_.store( row_, j             , matrix_.load(row_,j             ) + it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size    , matrix_.load(row_,j+IT::size    ) + it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*2UL, matrix_.load(row_,j+IT::size*2UL) + it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*3UL, matrix_.load(row_,j+IT::size*3UL) + it.load() ); it += IT::size;
   }
   for( ; j<jpos; j+=IT::size, it+=IT::size ) {
      matrix_.store( row_, j, matrix_.load(row_,j) + it.load() );
   }
   for( ; remainder && j<columns; ++j, ++it ) {
      matrix_(row_,j) += *it;
   }
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO,SF>::addAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) += element->value();
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseRow<MT,SO,SF>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) -= (~rhs)[j    ];
      matrix_(row_,j+1UL) -= (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) -= (~rhs)[jpos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseRow<MT,SO,SF>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t columns( size() );

   const size_t jpos( ( remainder )?( columns & size_t(-IT::size) ):( columns ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( columns - ( columns % (IT::size) ) ) == jpos, "Invalid end calculation" );

   size_t j( 0UL );
   typename VT::ConstIterator it( (~rhs).begin() );

   for( ; (j+IT::size*3UL) < jpos; j+=IT::size*4UL ) {
      matrix_.store( row_, j             , matrix_.load(row_,j             ) - it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size    , matrix_.load(row_,j+IT::size    ) - it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*2UL, matrix_.load(row_,j+IT::size*2UL) - it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*3UL, matrix_.load(row_,j+IT::size*3UL) - it.load() ); it += IT::size;
   }
   for( ; j<jpos; j+=IT::size, it+=IT::size ) {
      matrix_.store( row_, j, matrix_.load(row_,j) - it.load() );
   }
   for( ; remainder && j<columns; ++j, ++it ) {
      matrix_(row_,j) -= *it;
   }
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
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO,SF>::subAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) -= element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseRow<MT,SO,SF>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) *= (~rhs)[j    ];
      matrix_(row_,j+1UL) *= (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) *= (~rhs)[jpos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO,SF>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseRow<MT,SO,SF>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t columns( size() );

   const size_t jpos( ( remainder )?( columns & size_t(-IT::size) ):( columns ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( columns - ( columns % (IT::size) ) ) == jpos, "Invalid end calculation" );

   size_t j( 0UL );
   typename VT::ConstIterator it( (~rhs).begin() );

   for( ; (j+IT::size*3UL) < jpos; j+=IT::size*4UL ) {
      matrix_.store( row_, j             , matrix_.load(row_,j             ) * it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size    , matrix_.load(row_,j+IT::size    ) * it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*2UL, matrix_.load(row_,j+IT::size*2UL) * it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*3UL, matrix_.load(row_,j+IT::size*3UL) * it.load() ); it += IT::size;
   }
   for( ; j<jpos; j+=IT::size, it+=IT::size ) {
      matrix_.store( row_, j, matrix_.load(row_,j) * it.load() );
   }
   for( ; remainder && j<columns; ++j, ++it ) {
      matrix_(row_,j) *= *it;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF >      // Symmetry flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO,SF>::multAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) = tmp[element->index()] * element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseRow for general column-major matrices.
// \ingroup dense_row
//
// This specialization of DenseRow adapts the class template to the requirements of general
// column-major matrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseRow<MT,false,false> : public DenseVector< DenseRow<MT,false,false>, true >
                               , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseRow<MT,false,false>            This;           //!< Type of this DenseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseRow&                     CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant row value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant row value.
   typedef typename If< Or< IsConst<MT>, Not< HasMutableDataAccess<MT> > >
                      , ConstPointer, ElementType* >::Type  Pointer;
   //**********************************************************************************************

   //**RowIterator class definition****************************************************************
   /*!\brief Iterator over the elements of the dense row.
   */
   template< typename MatrixType >  // Type of the dense matrix
   class RowIterator
   {
    private:
      //*******************************************************************************************
      //! Compilation switch for the return type of the value member function.
      /*! The \a returnConst compile time constant expression represents a compilation switch for
          the return type of the member access operators. In case the given matrix type \a MatrixType
          is const qualified, \a returnConst will be set to 1 and the member access operators will
          return a reference to const. Otherwise \a returnConst will be set to 0 and the member
          access operators will offer write access to the dense matrix elements. */
      enum { returnConst = IsConst<MatrixType>::value };
      //*******************************************************************************************

    public:
      //**Type definitions*************************************************************************
      //! Return type for the access to the value of a dense element.
      typedef typename IfTrue< returnConst
                             , typename MatrixType::ConstReference
                             , typename MatrixType::Reference >::Type  Reference;

      typedef std::random_access_iterator_tag  IteratorCategory;  //!< The iterator category.
      typedef RemoveReference<Reference>       ValueType;         //!< Type of the underlying elements.
      typedef ValueType*                       PointerType;       //!< Pointer return type.
      typedef Reference                        ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                        DifferenceType;    //!< Difference between two iterators.

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
         : matrix_( NULL )  // The dense matrix containing the row.
         , row_   ( 0UL  )  // The current row index.
         , column_( 0UL  )  // The current column index.
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
         : matrix_( &matrix )  // The dense matrix containing the row.
         , row_   ( row     )  // The current row index.
         , column_( column  )  // The current column index.
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different RowIterator instances.
      //
      // \param it The row iterator to be copied.
      */
      template< typename MatrixType2 >
      inline RowIterator( const RowIterator<MatrixType2>& it )
         : matrix_( it.matrix_ )  // The dense matrix containing the row.
         , row_   ( it.row_    )  // The current row index.
         , column_( it.column_ )  // The current column index.
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline RowIterator& operator+=( size_t inc ) {
         column_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline RowIterator& operator-=( size_t dec ) {
         column_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline RowIterator& operator++() {
         ++column_;
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

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline RowIterator& operator--() {
         --column_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const RowIterator operator--( int ) {
         const RowIterator tmp( *this );
         --(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Subscript operator***********************************************************************
      /*!\brief Direct access to the dense row elements.
      //
      // \param index Access index.
      // \return Reference to the accessed value.
      */
      inline ReferenceType operator[]( size_t index ) const {
         return (*matrix_)(row_,column_+index);
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense row element at the current iterator position.
      //
      // \return Reference to the current value.
      */
      inline ReferenceType operator*() const {
         return (*matrix_)(row_,column_);
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense row element at the current iterator position.
      //
      // \return Pointer to the dense row element at the current iterator position.
      */
      inline PointerType operator->() const {
         return &(*matrix_)(row_,column_);
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2 >
      inline bool operator==( const RowIterator<MatrixType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ == rhs.column_ );
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2 >
      inline bool operator!=( const RowIterator<MatrixType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      template< typename MatrixType2 >
      inline bool operator<( const RowIterator<MatrixType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ < rhs.column_ );
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      template< typename MatrixType2 >
      inline bool operator>( const RowIterator<MatrixType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ > rhs.column_ );
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      template< typename MatrixType2 >
      inline bool operator<=( const RowIterator<MatrixType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ <= rhs.column_ );
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      template< typename MatrixType2 >
      inline bool operator>=( const RowIterator<MatrixType2>& rhs ) const {
         return ( matrix_ == rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ >= rhs.column_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two row iterators.
      //
      // \param rhs The right-hand side row iterator.
      // \return The number of elements between the two row iterators.
      */
      inline DifferenceType operator-( const RowIterator& rhs ) const {
         return column_ - rhs.column_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a RowIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const RowIterator operator+( const RowIterator& it, size_t inc ) {
         return RowIterator( *it.matrix_, it.row_, it.column_+inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a RowIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const RowIterator operator+( size_t inc, const RowIterator& it ) {
         return RowIterator( *it.matrix_, it.row_, it.column_+inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a RowIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const RowIterator operator-( const RowIterator& it, size_t dec ) {
         return RowIterator( *it.matrix_, it.row_, it.column_-dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType* matrix_;  //!< The dense matrix containing the row.
      size_t      row_;     //!< The current row index.
      size_t      column_;  //!< The current column index.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2 > friend class RowIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef RowIterator<const MT>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, RowIterator<MT> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseRow( MT& matrix, size_t index );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseRow& operator=( const ElementType& rhs );
   inline DenseRow& operator=( const DenseRow& rhs );

   template< typename VT > inline DenseRow& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator*=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline DenseRow& operator*=( const SparseVector<VT,true>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseRow >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseRow >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t    size() const;
                              inline size_t    capacity() const;
                              inline size_t    nonZeros() const;
                              inline void      reset();
   template< typename Other > inline DenseRow& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool SO2, bool SF2 >
   inline bool canAlias( const DenseRow<MT2,SO2,SF2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool SO2, bool SF2 >
   inline bool isAliased( const DenseRow<MT2,SO2,SF2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   template< typename VT > inline void assign    ( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void assign    ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void addAssign ( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void addAssign ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void subAssign ( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void subAssign ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline void multAssign( const DenseVector <VT,true>& rhs );
   template< typename VT > inline void multAssign( const SparseVector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The dense matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool SF2 > friend class DenseRow;

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const DenseRow<MT2,SO2,SF2>& row );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const DenseRow<MT2,SO2,SF2>& a, const DenseRow<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< DenseRow<MT2,SO2,SF2> >::Type
      derestrict( DenseRow<MT2,SO2,SF2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE        ( MT );
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
/*!\brief The constructor for DenseRow.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false,false>::DenseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the row
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::Reference
   DenseRow<MT,false,false>::operator[]( size_t index )
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstReference
   DenseRow<MT,false,false>::operator[]( size_t index ) const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::Reference
   DenseRow<MT,false,false>::at( size_t index )
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstReference
   DenseRow<MT,false,false>::at( size_t index ) const
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
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::Pointer DenseRow<MT,false,false>::data()
{
   return matrix_.data() + row_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstPointer DenseRow<MT,false,false>::data() const
{
   return matrix_.data() + row_;
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::Iterator DenseRow<MT,false,false>::begin()
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstIterator DenseRow<MT,false,false>::begin() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstIterator DenseRow<MT,false,false>::cbegin() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::Iterator DenseRow<MT,false,false>::end()
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstIterator DenseRow<MT,false,false>::end() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,false>::ConstIterator DenseRow<MT,false,false>::cend() const
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
/*!\brief Homogenous assignment to all row elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row.
//
// This function homogeneously assigns the given value to all elements of the row. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator=( const ElementType& rhs )
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

   for( size_t j=jbegin; j<jend; ++j )
      matrix_(row_,j) = rhs;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseRow.
//
// \param rhs Dense row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator=( const DenseRow& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, rhs );

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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

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
      const ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAddAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !trySubAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a dense vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator*=( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryMultAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::operator*=( const SparseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const ResultType right( *this * (~rhs) );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, right );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a dense row and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the vector.
//
// This operator cannot be used for rows on lower or upper unitriangular matrices. The attempt
// to scale such a row results in a compilation error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,false,false> >::Type&
   DenseRow<MT,false,false>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   return operator=( (*this) * rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense row by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the vector.
//
// This operator cannot be used for rows on lower or upper unitriangular matrices. The attempt
// to scale such a row results in a compilation error!
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,false,false> >::Type&
   DenseRow<MT,false,false>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   return operator=( (*this) / rhs );
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false,false>::size() const
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row.
//
// \return The capacity of the dense row.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false,false>::capacity() const
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false,false>::nonZeros() const
{
   const size_t columns( size() );
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns; ++j )
      if( !isDefault( matrix_(row_,j) ) )
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
inline void DenseRow<MT,false,false>::reset()
{
   using blaze::clear;

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

   for( size_t j=jbegin; j<jend; ++j )
      clear( matrix_(row_,j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseRow<MT,false,false>& DenseRow<MT,false,false>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jbegin( ( IsUpper<MT>::value )
                        ?( ( IsStrictlyUpper<MT>::value )
                           ?( row_+1UL )
                           :( row_ ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower<MT>::value )
                        ?( ( IsStrictlyLower<MT>::value )
                           ?( row_ )
                           :( row_+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_(row_,j) *= scalar;
   }

   return *this;
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
/*!\brief Returns whether the dense row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,false,false>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can alias with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense row
        , bool SO2       // Storage order of the foreign dense row
        , bool SF2 >     // Symmetry flag of the foreign dense row
inline bool DenseRow<MT,false,false>::canAlias( const DenseRow<MT2,SO2,SF2>* alias ) const
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row_ == alias->row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,false,false>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense row
        , bool SO2       // Storage order of the foreign dense row
        , bool SF2 >     // Symmetry flag of the foreign dense row
inline bool DenseRow<MT,false,false>::isAliased( const DenseRow<MT2,SO2,SF2>* alias ) const
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row_ == alias->row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is properly aligned in memory.
//
// \return \a true in case the dense row is aligned, \a false if not.
//
// This function returns whether the dense row is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense row are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseRow<MT,false,false>::isAligned() const
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can be used in SMP assignments.
//
// \return \a true in case the dense row can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense row).
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseRow<MT,false,false>::canSMPAssign() const
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void DenseRow<MT,false,false>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) = (~rhs)[j    ];
      matrix_(row_,j+1UL) = (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) = (~rhs)[jpos];
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,false>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) = element->value();
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void DenseRow<MT,false,false>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) += (~rhs)[j    ];
      matrix_(row_,j+1UL) += (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) += (~rhs)[jpos];
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,false>::addAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) += element->value();
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void DenseRow<MT,false,false>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) -= (~rhs)[j    ];
      matrix_(row_,j+1UL) -= (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) -= (~rhs)[jpos];
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,false>::subAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void DenseRow<MT,false,false>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jpos( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jpos; j+=2UL ) {
      matrix_(row_,j    ) *= (~rhs)[j    ];
      matrix_(row_,j+1UL) *= (~rhs)[j+1UL];
   }
   if( jpos < (~rhs).size() )
      matrix_(row_,jpos) *= (~rhs)[jpos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,false>::multAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) = tmp[element->index()] * element->value();
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
/*!\brief Specialization of DenseRow for symmetric column-major matrices.
// \ingroup dense_row
//
// This specialization of DenseRow adapts the class template to the requirements of symmetric
// column-major matrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseRow<MT,false,true> : public DenseVector< DenseRow<MT,false,true>, true >
                              , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename If< IsExpression<MT>, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the row element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseRow<MT,false,true>             This;           //!< Type of this DenseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseRow&                     CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
   typedef typename If< IsConst<MT>, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant row value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant row value.
   typedef typename If< Or< IsConst<MT>, Not< HasMutableDataAccess<MT> > >
                      , ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<MT>, ConstIterator, typename MT::Iterator >::Type  Iterator;
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
   explicit inline DenseRow( MT& matrix, size_t index );
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
   inline Pointer        data  ();
   inline ConstPointer   data  () const;
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
   inline DenseRow& operator=( const ElementType& rhs );
   inline DenseRow& operator=( const DenseRow& rhs );

   template< typename VT > inline DenseRow& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator*=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline DenseRow& operator*=( const SparseVector<VT,true>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseRow >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseRow >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t    size() const;
                              inline size_t    capacity() const;
                              inline size_t    nonZeros() const;
                              inline void      reset();
   template< typename Other > inline DenseRow& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedAddAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedSubAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedMultAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::multiplication };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename MT2, bool SO2, bool SF2 >
   inline bool canAlias( const DenseRow<MT2,SO2,SF2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename MT2, bool SO2, bool SF2 >
   inline bool isAliased( const DenseRow<MT2,SO2,SF2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   BLAZE_ALWAYS_INLINE IntrinsicType load ( size_t index ) const;
   BLAZE_ALWAYS_INLINE IntrinsicType loada( size_t index ) const;
   BLAZE_ALWAYS_INLINE IntrinsicType loadu( size_t index ) const;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void storea( size_t index, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const IntrinsicType& value );
   BLAZE_ALWAYS_INLINE void stream( size_t index, const IntrinsicType& value );

   template< typename VT >
   inline typename DisableIf< VectorizedAssign<VT> >::Type
      assign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedAssign<VT> >::Type
      assign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void assign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedAddAssign<VT> >::Type
      addAssign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedAddAssign<VT> >::Type
      addAssign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void addAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedSubAssign<VT> >::Type
      subAssign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedSubAssign<VT> >::Type
      subAssign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void subAssign( const SparseVector<VT,true>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedMultAssign<VT> >::Type
      multAssign( const DenseVector<VT,true>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedMultAssign<VT> >::Type
      multAssign( const DenseVector<VT,true>& rhs );

   template< typename VT > inline void multAssign( const SparseVector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The dense matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, bool SO2, bool SF2 > friend class DenseRow;

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isIntact( const DenseRow<MT2,SO2,SF2>& row );

   template< typename MT2, bool SO2, bool SF2 >
   friend bool isSame( const DenseRow<MT2,SO2,SF2>& a, const DenseRow<MT2,SO2,SF2>& b );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryAddAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool trySubAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2, typename VT >
   friend bool tryMultAssign( const DenseRow<MT2,SO2,SF2>& lhs, const Vector<VT,true>& rhs, size_t index );

   template< typename MT2, bool SO2, bool SF2 >
   friend typename DerestrictTrait< DenseRow<MT2,SO2,SF2> >::Type
      derestrict( DenseRow<MT2,SO2,SF2>& dm );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
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
/*!\brief The constructor for DenseRow.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false,true>::DenseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the row
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::Reference
   DenseRow<MT,false,true>::operator[]( size_t index )
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstReference
   DenseRow<MT,false,true>::operator[]( size_t index ) const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::Reference
   DenseRow<MT,false,true>::at( size_t index )
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstReference
   DenseRow<MT,false,true>::at( size_t index ) const
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
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::Pointer DenseRow<MT,false,true>::data()
{
   return matrix_.data( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row. Note that in case
// of a column-major matrix you can NOT assume that the row elements lie adjacent to each other!
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstPointer DenseRow<MT,false,true>::data() const
{
   return matrix_.data( row_ );
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::Iterator DenseRow<MT,false,true>::begin()
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstIterator DenseRow<MT,false,true>::begin() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstIterator DenseRow<MT,false,true>::cbegin() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::Iterator DenseRow<MT,false,true>::end()
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstIterator DenseRow<MT,false,true>::end() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false,true>::ConstIterator DenseRow<MT,false,true>::cend() const
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
/*!\brief Homogenous assignment to all row elements.
//
// \param rhs Scalar value to be assigned to all row elements.
// \return Reference to the assigned row.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator=( const ElementType& rhs )
{
   const size_t ibegin( ( IsLower<MT>::value )
                        ?( ( IsUniLower<MT>::value || IsStrictlyLower<MT>::value )
                           ?( row_+1UL )
                           :( row_ ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper<MT>::value )
                        ?( ( IsUniUpper<MT>::value || IsStrictlyUpper<MT>::value )
                           ?( row_ )
                           :( row_+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i )
      matrix_(i,row_) = rhs;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseRow.
//
// \param rhs Dense row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator=( const DenseRow& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, rhs );

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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator=( const Vector<VT,true>& rhs )
{
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
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryAddAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !trySubAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a dense vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator*=( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<MT>, typename VT::CompositeType, const VT& >::Type  Right;
   Right right( ~rhs );

   if( !tryMultAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::operator*=( const SparseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const ResultType right( *this * (~rhs) );

   if( !tryAssign( matrix_, right, row_, 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, right );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a dense row and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the vector.
//
// This operator cannot be used for rows on lower or upper unitriangular matrices. The attempt
// to scale such a row results in a compilation error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,false,true> >::Type&
   DenseRow<MT,false,true>::operator*=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   return operator=( (*this) * rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense row by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the vector.
//
// This operator cannot be used for rows on lower or upper unitriangular matrices. The attempt
// to scale such a row results in a compilation error!
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,false,true> >::Type&
   DenseRow<MT,false,true>::operator/=( Other rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   return operator=( (*this) / rhs );
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false,true>::size() const
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense row.
//
// \return The capacity of the dense row.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false,true>::capacity() const
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false,true>::nonZeros() const
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
template< typename MT >  // Type of the dense matrix
inline void DenseRow<MT,false,true>::reset()
{
   matrix_.reset( row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
//
// This function scales all elements of the row by the given scalar value \a scalar. Note that
// the function cannot be used to scale a row on a lower or upper unitriangular matrix. The
// attempt to scale such a row results in a compile time error!
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseRow<MT,false,true>& DenseRow<MT,false,true>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t ibegin( ( IsLower<MT>::value )
                        ?( ( IsStrictlyLower<MT>::value )
                           ?( row_+1UL )
                           :( row_ ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper<MT>::value )
                        ?( ( IsStrictlyUpper<MT>::value )
                           ?( row_ )
                           :( row_+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_(i,row_) *= scalar;
   }

   return *this;
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
/*!\brief Returns whether the dense row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,false,true>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can alias with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address can alias with the dense row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense row
        , bool SO2       // Storage order of the foreign dense row
        , bool SF2 >     // Symmetry flag of the foreign dense row
inline bool DenseRow<MT,false,true>::canAlias( const DenseRow<MT2,SO2,SF2>* alias ) const
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row_ == alias->row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,false,true>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is aliased with the given dense row \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense row, \a false if not.
//
// This function returns whether the given address is aliased with the dense row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >  // Type of the dense matrix
template< typename MT2   // Data type of the foreign dense row
        , bool SO2       // Storage order of the foreign dense row
        , bool SF2 >     // Symmetry flag of the foreign dense row
inline bool DenseRow<MT,false,true>::isAliased( const DenseRow<MT2,SO2,SF2>* alias ) const
{
   return matrix_.isAliased( &alias->matrix_ ) && ( row_ == alias->row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row is properly aligned in memory.
//
// \return \a true in case the dense row is aligned, \a false if not.
//
// This function returns whether the dense row is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the dense row are guaranteed to conform to the
// alignment restrictions of the element type \a Type.
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseRow<MT,false,true>::isAligned() const
{
   return matrix_.isAligned();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense row can be used in SMP assignments.
//
// \return \a true in case the dense row can be used in SMP assignments, \a false if not.
//
// This function returns whether the dense row can be used in SMP assignments. In contrast to
// the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size
// of the dense row).
*/
template< typename MT >  // Type of the dense matrix
inline bool DenseRow<MT,false,true>::canSMPAssign() const
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded intrinsic element.
//
// This function performs a load of a specific intrinsic element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseRow<MT,false,true>::IntrinsicType
   DenseRow<MT,false,true>::load( size_t index ) const
{
   return matrix_.load( index, row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseRow<MT,false,true>::IntrinsicType
   DenseRow<MT,false,true>::loada( size_t index ) const
{
   return matrix_.loada( index, row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE typename DenseRow<MT,false,true>::IntrinsicType
   DenseRow<MT,false,true>::loadu( size_t index ) const
{
   return matrix_.loadu( index, row_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs a store a specific intrinsic element of the dense row. The index
// must be smaller than the number of matrix columns. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void DenseRow<MT,false,true>::store( size_t index, const IntrinsicType& value )
{
   matrix_.store( index, row_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void DenseRow<MT,false,true>::storea( size_t index, const IntrinsicType& value )
{
   matrix_.storea( index, row_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unligned store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store a specific intrinsic element of the dense row.
// The index must be smaller than the number of matrix columns. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void DenseRow<MT,false,true>::storeu( size_t index, const IntrinsicType& value )
{
   matrix_.storeu( index, row_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of an intrinsic element of the dense row.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific intrinsic element of the
// dense row. The index must be smaller than the number of matrix columns. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized evaluation
// of expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT >  // Type of the dense matrix
BLAZE_ALWAYS_INLINE void DenseRow<MT,false,true>::stream( size_t index, const IntrinsicType& value )
{
   matrix_.stream( index, row_, value );
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseRow<MT,false,true>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row_) = (~rhs)[i    ];
      matrix_(i+1UL,row_) = (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      matrix_(ipos,row_) = (~rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseRow<MT,false,true>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t rows( size() );

   const size_t ipos( ( remainder )?( rows & size_t(-IT::size) ):( rows ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( rows - ( rows % (IT::size) ) ) == ipos, "Invalid end calculation" );

   if( useStreaming && rows > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(~rhs).isAliased( &matrix_ ) )
   {
      size_t i( 0UL );

      for( ; i<ipos; i+=IT::size ) {
         matrix_.stream( i, row_, (~rhs).load(i) );
      }
      for( ; remainder && i<rows; ++i ) {
         matrix_(i,row_) = (~rhs)[i];
      }
   }
   else
   {
      size_t i( 0UL );
      typename VT::ConstIterator it( (~rhs).begin() );

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
         matrix_.store( i             , row_, it.load() ); it += IT::size;
         matrix_.store( i+IT::size    , row_, it.load() ); it += IT::size;
         matrix_.store( i+IT::size*2UL, row_, it.load() ); it += IT::size;
         matrix_.store( i+IT::size*3UL, row_, it.load() ); it += IT::size;
      }
      for( ; i<ipos; i+=IT::size, it+=IT::size ) {
         matrix_.store( i, row_, it.load() );
      }
      for( ; remainder && i<rows; ++i, ++it ) {
         matrix_(i,row_) = *it;
      }
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,true>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),row_) = element->value();
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseRow<MT,false,true>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row_) += (~rhs)[i    ];
      matrix_(i+1UL,row_) += (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      matrix_(ipos,row_) += (~rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseRow<MT,false,true>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t rows( size() );

   const size_t ipos( ( remainder )?( rows & size_t(-IT::size) ):( rows ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( rows - ( rows % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      matrix_.store( i             , row_, matrix_.load(i             ,row_) + it.load() ); it += IT::size;
      matrix_.store( i+IT::size    , row_, matrix_.load(i+IT::size    ,row_) + it.load() ); it += IT::size;
      matrix_.store( i+IT::size*2UL, row_, matrix_.load(i+IT::size*2UL,row_) + it.load() ); it += IT::size;
      matrix_.store( i+IT::size*3UL, row_, matrix_.load(i+IT::size*3UL,row_) + it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      matrix_.store( i, row_, matrix_.load(i,row_) + it.load() );
   }
   for( ; remainder && i<rows; ++i, ++it ) {
      matrix_(i,row_) += *it;
   }
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,true>::addAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),row_) += element->value();
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseRow<MT,false,true>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row_) -= (~rhs)[i    ];
      matrix_(i+1UL,row_) -= (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      matrix_(ipos,row_) -= (~rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseRow<MT,false,true>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t rows( size() );

   const size_t ipos( ( remainder )?( rows & size_t(-IT::size) ):( rows ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( rows - ( rows % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      matrix_.store( i             , row_, matrix_.load(i             ,row_) - it.load() ); it += IT::size;
      matrix_.store( i+IT::size    , row_, matrix_.load(i+IT::size    ,row_) - it.load() ); it += IT::size;
      matrix_.store( i+IT::size*2UL, row_, matrix_.load(i+IT::size*2UL,row_) - it.load() ); it += IT::size;
      matrix_.store( i+IT::size*3UL, row_, matrix_.load(i+IT::size*3UL,row_) - it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      matrix_.store( i, row_, matrix_.load(i,row_) - it.load() );
   }
   for( ; remainder && i<rows; ++i, ++it ) {
      matrix_(i,row_) -= *it;
   }
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
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,true>::subAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),row_) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseRow<MT,false,true>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      matrix_(i    ,row_) *= (~rhs)[i    ];
      matrix_(i+1UL,row_) *= (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      matrix_(ipos,row_) *= (~rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,false,true>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseRow<MT,false,true>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const bool remainder( !IsPadded<MT>::value || !IsPadded<VT>::value );
   const size_t rows( size() );

   const size_t ipos( ( remainder )?( rows & size_t(-IT::size) ):( rows ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( rows - ( rows % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      matrix_.store( i             , row_, matrix_.load(i             ,row_) * it.load() ); it += IT::size;
      matrix_.store( i+IT::size    , row_, matrix_.load(i+IT::size    ,row_) * it.load() ); it += IT::size;
      matrix_.store( i+IT::size*2UL, row_, matrix_.load(i+IT::size*2UL,row_) * it.load() ); it += IT::size;
      matrix_.store( i+IT::size*3UL, row_, matrix_.load(i+IT::size*3UL,row_) * it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      matrix_.store( i, row_, matrix_.load(i,row_) * it.load() );
   }
   for( ; remainder && i<rows; ++i, ++it ) {
      matrix_(i,row_) *= *it;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,false,true>::multAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),row_) = tmp[element->index()] * element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  DENSEROW OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseRow operators */
//@{
template< typename MT, bool SO, bool SF >
inline void reset( DenseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline void clear( DenseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline bool isDefault( const DenseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline bool isIntact( const DenseRow<MT,SO,SF>& row );

template< typename MT, bool SO, bool SF >
inline bool isSame( const DenseRow<MT,SO,SF>& a, const DenseRow<MT,SO,SF>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dense row.
// \ingroup dense_row
//
// \param row The dense row to be resetted.
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void reset( DenseRow<MT,SO,SF>& row )
{
   row.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dense row.
// \ingroup dense_row
//
// \param row The dense row to be cleared.
// \return void
//
// Clearing a dense row is equivalent to resetting it via the reset() function.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline void clear( DenseRow<MT,SO,SF>& row )
{
   row.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense row is in default state.
// \ingroup dense_row
//
// \param row The dense row to be tested for its default state.
// \return \a true in case the given dense row is component-wise zero, \a false otherwise.
//
// This function checks whether the dense row is in default state. For instance, in case the
// row is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all row elements are 0 and \a false in case any row element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isDefault( const DenseRow<MT,SO,SF>& row )
{
   for( size_t i=0UL; i<row.size(); ++i )
      if( !isDefault( row[i] ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given dense row are intact.
// \ingroup dense_row
//
// \param row The dense row to be tested.
// \return \a true in case the given row's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dense row are intact, i.e. if its state
// is valid. In case the invariants are intact, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isIntact( const DenseRow<MT,SO,SF>& row )
{
   return ( row.row_ <= row.matrix_.rows() &&
            isIntact( row.matrix_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given dense rows represent the same observable state.
// \ingroup dense_row
//
// \param a The first dense row to be tested for its state.
// \param b The second dense row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given dense rows refer to exactly the
// same range of the same dense matrix. In case both rows represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline bool isSame( const DenseRow<MT,SO,SF>& a, const DenseRow<MT,SO,SF>& b )
{
   return ( isSame( a.matrix_, b.matrix_ ) && ( a.row_ == b.row_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a dense row.
// \ingroup dense_row
//
// \param lhs The target left-hand side dense row.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool tryAssign( const DenseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a dense row.
// \ingroup dense_row
//
// \param lhs The target left-hand side dense row.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool tryAddAssign( const DenseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a dense row.
// \ingroup dense_row
//
// \param lhs The target left-hand side dense row.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool trySubAssign( const DenseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a dense row.
// \ingroup dense_row
//
// \param lhs The target left-hand side dense row.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , bool SF        // Symmetry flag
        , typename VT >  // Type of the right-hand side vector
inline bool tryMultAssign( const DenseRow<MT,SO,SF>& lhs, const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.matrix_, ~rhs, lhs.row_, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given dense row.
// \ingroup dense_row
//
// \param row The dense row to be derestricted.
// \return Dense row without access restrictions.
//
// This function removes all restrictions on the data access to the given dense row. It returns a
// row object that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , bool SF >    // Symmetry flag
inline typename DerestrictTrait< DenseRow<MT,SO,SF> >::Type
   derestrict( DenseRow<MT,SO,SF>& row )
{
   typedef typename DerestrictTrait< DenseRow<MT,SO,SF> >::Type  ReturnType;
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
struct IsRestricted< DenseRow<MT,SO,SF> > : public IsTrue< IsRestricted<MT>::value >
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
struct DerestrictTrait< DenseRow<MT,SO,SF> >
{
   typedef DenseRow< typename RemoveReference< typename DerestrictTrait<MT>::Type >::Type >  Type;
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
template< typename MT, bool SO, bool SF >
struct HasConstDataAccess< DenseRow<MT,SO,SF> >
   : public IsTrue< HasConstDataAccess<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF >
struct HasMutableDataAccess< DenseRow<MT,SO,SF> >
   : public IsTrue< HasMutableDataAccess<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF >
struct IsAligned< DenseRow<MT,SO,SF> >
   : public IsTrue< And< IsAligned<MT>, Or< IsRowMajorMatrix<MT>, IsSymmetric<MT> > >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF >
struct IsPadded< DenseRow<MT,SO,SF> >
   : public IsTrue< And< IsPadded<MT>, Or< IsRowMajorMatrix<MT>, IsSymmetric<MT> > >::value >
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
template< typename MT, bool SO, bool SF, typename T >
struct AddTrait< DenseRow<MT,SO,SF>, T >
{
   typedef typename AddTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct AddTrait< T, DenseRow<MT,SO,SF> >
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
struct SubTrait< DenseRow<MT,SO,SF>, T >
{
   typedef typename SubTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct SubTrait< T, DenseRow<MT,SO,SF> >
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
struct MultTrait< DenseRow<MT,SO,SF>, T >
{
   typedef typename MultTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct MultTrait< T, DenseRow<MT,SO,SF> >
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
struct CrossTrait< DenseRow<MT,SO,SF>, T >
{
   typedef typename CrossTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct CrossTrait< T, DenseRow<MT,SO,SF> >
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
struct DivTrait< DenseRow<MT,SO,SF>, T >
{
   typedef typename DivTrait< typename RowTrait<MT>::Type, T >::Type  Type;
};

template< typename T, typename MT, bool SO, bool SF >
struct DivTrait< T, DenseRow<MT,SO,SF> >
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
struct SubvectorTrait< DenseRow<MT,SO,SF> >
{
   typedef typename SubvectorTrait< typename DenseRow<MT,SO,SF>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
