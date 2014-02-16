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
#include <stdexcept>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Row.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Streaming.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>


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
   template< typename MT, bool SO >
   class DenseRow;
   \endcode

//  - MT: specifies the type of the dense matrix primitive. DenseRow can be used with every dense
//        matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the dense matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
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
template< typename MT                              // Type of the dense matrix
        , bool SO = IsRowMajorMatrix<MT>::value >  // Storage order
class DenseRow : public DenseVector< DenseRow<MT,SO>, true >
               , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, MT, MT& >::Type  Operand;

   //! Intrinsic trait for the row element type.
   typedef IntrinsicTrait<typename MT::ElementType>  IT;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense row will return
       references and iterators to const. Otherwise \a useConst will be set to 0 and the
       dense row will offer write access to the dense matrix elements both via the subscript
       operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseRow<MT,SO>                     This;           //!< Type of this DenseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseRow&                     CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Pointer to a constant row value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a constant row value.
   typedef typename SelectType< useConst, ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, typename MT::Iterator >::Type  Iterator;
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
                           inline DenseRow& operator= ( const ElementType& rhs );
                           inline DenseRow& operator= ( const DenseRow& rhs );
   template< typename VT > inline DenseRow& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator*=( const Vector<VT,true>& rhs );

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
      enum { value = vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct VectorizedAddAssign {
      enum { value = vectorizable && VT::vectorizable &&
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
      enum { value = vectorizable && VT::vectorizable &&
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
      enum { value = vectorizable && VT::vectorizable &&
                     IsSame<ElementType,typename VT::ElementType>::value &&
                     IntrinsicTrait<ElementType>::multiplication };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load  ( size_t index ) const;
   inline IntrinsicType loadu ( size_t index ) const;
   inline void          store ( size_t index, const IntrinsicType& value );
   inline void          storeu( size_t index, const IntrinsicType& value );
   inline void          stream( size_t index, const IntrinsicType& value );

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

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
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
        , bool SO >    // Storage order
inline DenseRow<MT,SO>::DenseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the row
   , row_   ( index  )  // The index of the row in the matrix
{
   if( matrix_.rows() <= index )
      throw std::invalid_argument( "Invalid row access index" );
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
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::Reference DenseRow<MT,SO>::operator[]( size_t index )
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
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::ConstReference DenseRow<MT,SO>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row_,index);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::Pointer DenseRow<MT,SO>::data()
{
   return matrix_.data( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the row elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense row.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::ConstPointer DenseRow<MT,SO>::data() const
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::Iterator DenseRow<MT,SO>::begin()
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::ConstIterator DenseRow<MT,SO>::begin() const
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::ConstIterator DenseRow<MT,SO>::cbegin() const
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::Iterator DenseRow<MT,SO>::end()
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::ConstIterator DenseRow<MT,SO>::end() const
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::ConstIterator DenseRow<MT,SO>::cend() const
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
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseRow<MT,SO>& DenseRow<MT,SO>::operator=( const ElementType& rhs )
{
   const size_t columns( size() );

   for( size_t j=0UL; j<columns; ++j )
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
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseRow<MT,SO>& DenseRow<MT,SO>::operator=( const DenseRow& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Row sizes do not match" );

   const size_t columns( size() );

   for( size_t j=0UL; j<columns; ++j )
      matrix_(row_,j) = rhs[j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO>& DenseRow<MT,SO>::operator=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      smpAssign( *this, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      smpAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO>& DenseRow<MT,SO>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO>& DenseRow<MT,SO>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,SO>& DenseRow<MT,SO>::operator*=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) || IsSparseVector<VT>::value ) {
      const typename VT::ResultType tmp( ~rhs );
      smpMultAssign( *this, tmp );
   }
   else {
      smpMultAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a vector and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the vector.
*/
template< typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,SO> >::Type&
   DenseRow<MT,SO>::operator*=( Other rhs )
{
   return operator=( (*this) * rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a vector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the vector.
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,SO> >::Type&
   DenseRow<MT,SO>::operator/=( Other rhs )
{
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
        , bool SO >    // Storage order
inline size_t DenseRow<MT,SO>::size() const
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
        , bool SO >    // Storage order
inline size_t DenseRow<MT,SO>::capacity() const
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
        , bool SO >    // Storage order
inline size_t DenseRow<MT,SO>::nonZeros() const
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
        , bool SO >    // Storage order
inline void DenseRow<MT,SO>::reset()
{
   matrix_.reset( row_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
*/
template< typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline DenseRow<MT,SO>& DenseRow<MT,SO>::scale( const Other& scalar )
{
   for( size_t j=0UL; j<size(); ++j ) {
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,SO>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DenseRow<MT,SO>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
        , bool SO >    // Storage order
inline bool DenseRow<MT,SO>::isAligned() const
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
        , bool SO >    // Storage order
inline bool DenseRow<MT,SO>::canSMPAssign() const
{
   return ( size() > OPENMP_DVECASSIGN_THRESHOLD );
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::IntrinsicType DenseRow<MT,SO>::load( size_t index ) const
{
   return matrix_.load( row_, index );
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
        , bool SO >    // Storage order
inline typename DenseRow<MT,SO>::IntrinsicType DenseRow<MT,SO>::loadu( size_t index ) const
{
   return matrix_.loadu( row_, index );
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
        , bool SO >    // Storage order
inline void DenseRow<MT,SO>::store( size_t index, const IntrinsicType& value )
{
   matrix_.store( row_, index, value );
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
        , bool SO >    // Storage order
inline void DenseRow<MT,SO>::storeu( size_t index, const IntrinsicType& value )
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
        , bool SO >    // Storage order
inline void DenseRow<MT,SO>::stream( size_t index, const IntrinsicType& value )
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseRow<MT,SO>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) = (~rhs)[j    ];
      matrix_(row_,j+1UL) = (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) = (~rhs)[jend];
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseRow<MT,SO>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t columns( size() );

   if( useStreaming && columns > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(~rhs).isAliased( &matrix_ ) )
   {
      for( size_t j=0UL; j<columns; j+=IT::size ) {
         matrix_.stream( row_, j, (~rhs).load(j) );
      }
   }
   else
   {
      const size_t jend( columns & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( columns - ( columns % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

      typename VT::ConstIterator it( (~rhs).begin() );
      for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
         matrix_.store( row_, j             , it.load() ); it += IT::size;
         matrix_.store( row_, j+IT::size    , it.load() ); it += IT::size;
         matrix_.store( row_, j+IT::size*2UL, it.load() ); it += IT::size;
         matrix_.store( row_, j+IT::size*3UL, it.load() ); it += IT::size;
      }
      for( size_t j=jend; j<columns; j+=IT::size, it+=IT::size ) {
         matrix_.store( row_, j, it.load() );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO>::assign( const SparseVector<VT,true>& rhs )
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseRow<MT,SO>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) += (~rhs)[j    ];
      matrix_(row_,j+1UL) += (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) += (~rhs)[jend];
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseRow<MT,SO>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t columns( size() );

   const size_t jend( columns & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( columns - ( columns % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

   typename VT::ConstIterator it( (~rhs).begin() );
   for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
      matrix_.store( row_, j             , matrix_.load(row_,j             ) + it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size    , matrix_.load(row_,j+IT::size    ) + it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*2UL, matrix_.load(row_,j+IT::size*2UL) + it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*3UL, matrix_.load(row_,j+IT::size*3UL) + it.load() ); it += IT::size;
   }
   for( size_t j=jend; j<columns; j+=IT::size, it+=IT::size ) {
      matrix_.store( row_, j, matrix_.load(row_,j) + it.load() );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO>::addAssign( const SparseVector<VT,true>& rhs )
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseRow<MT,SO>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) -= (~rhs)[j    ];
      matrix_(row_,j+1UL) -= (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) -= (~rhs)[jend];
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseRow<MT,SO>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t columns( size() );

   const size_t jend( columns & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( columns - ( columns % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

   typename VT::ConstIterator it( (~rhs).begin() );
   for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
      matrix_.store( row_, j             , matrix_.load(row_,j             ) - it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size    , matrix_.load(row_,j+IT::size    ) - it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*2UL, matrix_.load(row_,j+IT::size*2UL) - it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*3UL, matrix_.load(row_,j+IT::size*3UL) - it.load() ); it += IT::size;
   }
   for( size_t j=jend; j<columns; j+=IT::size, it+=IT::size ) {
      matrix_.store( row_, j, matrix_.load(row_,j) - it.load() );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO>::subAssign( const SparseVector<VT,true>& rhs )
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseRow<MT,SO>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) *= (~rhs)[j    ];
      matrix_(row_,j+1UL) *= (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) *= (~rhs)[jend];
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseRow<MT,SO>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseRow<MT,SO>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t columns( size() );

   const size_t jend( columns & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( columns - ( columns % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

   typename VT::ConstIterator it( (~rhs).begin() );
   for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
      matrix_.store( row_, j             , matrix_.load(row_,j             ) * it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size    , matrix_.load(row_,j+IT::size    ) * it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*2UL, matrix_.load(row_,j+IT::size*2UL) * it.load() ); it += IT::size;
      matrix_.store( row_, j+IT::size*3UL, matrix_.load(row_,j+IT::size*3UL) * it.load() ); it += IT::size;
   }
   for( size_t j=jend; j<columns; j+=IT::size, it+=IT::size ) {
      matrix_.store( row_, j, matrix_.load(row_,j) * it.load() );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void DenseRow<MT,SO>::multAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( *this );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) = tmp[element->index()] * element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseRow for column-major matrices.
// \ingroup dense_row
//
// This specialization of DenseRow adapts the class template to the requirements of
// column-major matrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseRow<MT,false> : public DenseVector< DenseRow<MT,false>, true >
                         , private Row
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, MT, MT& >::Type  Operand;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense row will return
       references and iterators to const. Otherwise \a useConst will be set to 0 and the
       dense row will offer write access to the dense matrix elements both via the subscript
       operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseRow<MT,false>                  This;           //!< Type of this DenseRow instance.
   typedef typename RowTrait<MT>::Type         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the row elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const ResultType                    CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant row value.
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;
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
      typedef typename SelectType< returnConst, typename MT::ConstReference, typename MT::Reference >::Type  Reference;

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
      /*!\brief Constructor for the RowIterator class.
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
   typedef typename SelectType< useConst, ConstIterator, RowIterator<MT> >::Type  Iterator;
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
                           inline DenseRow& operator= ( const ElementType& rhs );
                           inline DenseRow& operator= ( const DenseRow& rhs );
   template< typename VT > inline DenseRow& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline DenseRow& operator*=( const Vector<VT,true>& rhs );

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
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

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

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
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
inline DenseRow<MT,false>::DenseRow( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the row
   , row_   ( index  )  // The index of the row in the matrix
{
   if( matrix_.rows() <= index )
      throw std::invalid_argument( "Invalid row access index" );
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
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false>::Reference DenseRow<MT,false>::operator[]( size_t index )
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
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseRow<MT,false>::ConstReference DenseRow<MT,false>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row_,index);
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
inline typename DenseRow<MT,false>::Iterator DenseRow<MT,false>::begin()
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
inline typename DenseRow<MT,false>::ConstIterator DenseRow<MT,false>::begin() const
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
inline typename DenseRow<MT,false>::ConstIterator DenseRow<MT,false>::cbegin() const
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
inline typename DenseRow<MT,false>::Iterator DenseRow<MT,false>::end()
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
inline typename DenseRow<MT,false>::ConstIterator DenseRow<MT,false>::end() const
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
inline typename DenseRow<MT,false>::ConstIterator DenseRow<MT,false>::cend() const
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
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false>& DenseRow<MT,false>::operator=( const ElementType& rhs )
{
   const size_t columns( size() );

   for( size_t j=0UL; j<columns; ++j )
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
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseRow<MT,false>& DenseRow<MT,false>::operator=( const DenseRow& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Row sizes do not match" );

   const size_t columns( size() );

   for( size_t j=0UL; j<columns; ++j )
      matrix_(row_,j) = rhs[j];

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false>& DenseRow<MT,false>::operator=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const ResultType tmp( ~rhs );
      smpAssign( *this, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      smpAssign( *this, ~rhs );
   }

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false>& DenseRow<MT,false>::operator+=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, ~rhs );
   }

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false>& DenseRow<MT,false>::operator-=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense row.
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseRow<MT,false>& DenseRow<MT,false>::operator*=( const Vector<VT,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) || IsSparseVector<VT>::value ) {
      const typename VT::ResultType tmp( ~rhs );
      smpMultAssign( *this, tmp );
   }
   else {
      smpMultAssign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a vector and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the vector.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,false> >::Type&
   DenseRow<MT,false>::operator*=( Other rhs )
{
   for( size_t j=0UL; j<size(); ++j )
      matrix_(row_,j) *= rhs;
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a vector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the vector.
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseRow<MT,false> >::Type&
   DenseRow<MT,false>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( size_t j=0UL; j<size(); ++j )
         matrix_(row_,j) *= tmp;
   }
   else {
      for( size_t j=0UL; j<size(); ++j )
         matrix_(row_,j) /= rhs;
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
template< typename MT >  // Type of the dense matrix
inline size_t DenseRow<MT,false>::size() const
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
inline size_t DenseRow<MT,false>::capacity() const
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
inline size_t DenseRow<MT,false>::nonZeros() const
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
inline void DenseRow<MT,false>::reset()
{
   using blaze::reset;
   const size_t columns( size() );
   for( size_t j=0UL; j<columns; ++j )
      reset( matrix_(row_,j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the dense row.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseRow<MT,false>& DenseRow<MT,false>::scale( const Other& scalar )
{
   for( size_t j=0UL; j<size(); ++j ) {
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
inline bool DenseRow<MT,false>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
inline bool DenseRow<MT,false>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
inline bool DenseRow<MT,false>::isAligned() const
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
inline bool DenseRow<MT,false>::canSMPAssign() const
{
   return ( size() > OPENMP_DVECASSIGN_THRESHOLD );
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
inline void DenseRow<MT,false>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) = (~rhs)[j    ];
      matrix_(row_,j+1UL) = (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) = (~rhs)[jend];
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
inline void DenseRow<MT,false>::assign( const SparseVector<VT,true>& rhs )
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
inline void DenseRow<MT,false>::addAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) += (~rhs)[j    ];
      matrix_(row_,j+1UL) += (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) += (~rhs)[jend];
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
inline void DenseRow<MT,false>::addAssign( const SparseVector<VT,true>& rhs )
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
inline void DenseRow<MT,false>::subAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) -= (~rhs)[j    ];
      matrix_(row_,j+1UL) -= (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) -= (~rhs)[jend];
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
inline void DenseRow<MT,false>::subAssign( const SparseVector<VT,true>& rhs )
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
inline void DenseRow<MT,false>::multAssign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t jend( (~rhs).size() & size_t(-2) );
   for( size_t j=0UL; j<jend; j+=2UL ) {
      matrix_(row_,j    ) *= (~rhs)[j    ];
      matrix_(row_,j+1UL) *= (~rhs)[j+1UL];
   }
   if( jend < (~rhs).size() )
      matrix_(row_,jend) *= (~rhs)[jend];
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
inline void DenseRow<MT,false>::multAssign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( *this );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(row_,element->index()) = tmp[element->index()] * element->value();
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
template< typename MT, bool SO >
inline void reset( DenseRow<MT,SO>& row );

template< typename MT, bool SO >
inline void clear( DenseRow<MT,SO>& row );

template< typename MT, bool SO >
inline bool isDefault( const DenseRow<MT,SO>& row );
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
        , bool SO >    // Storage order
inline void reset( DenseRow<MT,SO>& row )
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
        , bool SO >    // Storage order
inline void clear( DenseRow<MT,SO>& row )
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
        , bool SO >    // Storage order
inline bool isDefault( const DenseRow<MT,SO>& row )
{
   for( size_t i=0UL; i<row.size(); ++i )
      if( !isDefault( row[i] ) ) return false;
   return true;
}
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct SubvectorTrait< DenseRow<MT,SO> >
{
   typedef typename SubvectorTrait< typename DenseRow<MT,SO>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
