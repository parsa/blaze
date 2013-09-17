//=================================================================================================
/*!
//  \file blaze/math/views/SparseSubmatrix.h
//  \brief Header file for the SparseSubmatrix class template
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_VIEWS_SPARSESUBMATRIX_H_
#define _BLAZE_MATH_VIEWS_SPARSESUBMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <stdexcept>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsMatAbsExpr.h>
#include <blaze/math/typetraits/IsMatEvalExpr.h>
#include <blaze/math/typetraits/IsMatMatAddExpr.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsMatMatSubExpr.h>
#include <blaze/math/typetraits/IsMatScalarDivExpr.h>
#include <blaze/math/typetraits/IsMatScalarMultExpr.h>
#include <blaze/math/typetraits/IsMatTransExpr.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsVecTVecMultExpr.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
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
/*!\brief View to a specific submatrix of a sparse matrix.
// \ingroup sparse_submatrix
//
// The SparseSubmatrix template represents a view to a specific submatrix of a sparse matrix
// primitive. The type of the sparse matrix is specified via the first template parameter:

   \code
   template< typename MT, bool SO >
   class SparseSubmatrix;
   \endcode

//  - MT: specifies the type of the sparse matrix primitive. SparseSubmatrix can be used with any
//        sparse matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the sparse matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//
//
// \n \section sparse_submatrix_setup Setup of Sparse Submatrices
//
// A view to a sparse submatrix can very conveniently be created via the \c sub() function. This
// view can be treated as any other sparse matrix, i.e. it can be assigned to, it can be copied
// from, and it can be used in arithmetic operations. The view can also be used on both sides of
// an assignment: The submatrix can be either used as an alias to grant write access to a specific
// submatrix of a sparse matrix primitive on the left-hand side of an assignment or to grant
// read-access to a specific submatrix of a sparse matrix primitive or expression on the right-hand
// side of an assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;
   typedef blaze::DynamicMatrix<double,blaze::columnMajor>  DenseMatrixType;

   SparseMatrixType A, B;
   DenseMatrixType C;
   // ... Resizing and initialization

   // Creating a sparse submatrix of size 8x4, starting in row 0 and column 2
   blaze::SparseSubmatrix<SparseMatrixType> sm = sub( A, 0UL, 2UL, 8UL, 4UL );

   // Setting the submatrix of A to a 8x4 submatrix of B
   sm = sub( B, 0UL, 0UL, 8UL, 4UL );

   // Copying the dense matrix C into another 8x4 submatrix of A
   sub( A, 8UL, 2UL, 8UL, 4UL ) = C;

   // Assigning part of the result of a matrix addition to the first submatrix
   sm = sub( sub( B, 0UL, 0UL, 8UL, 4UL ) + C );
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
   blaze::SparseSubmatrix<MatrixType> sm = sub( A, 4UL, 4UL, 8UL, 8UL );

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
   SubmatrixType sm = sub( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( SubmatrixType::Iterator it=sm.begin(0); it!=sv.end(0); ++it ) {
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
   SubmatrixType sm( sub( A, 10UL, 10UL, 16UL, 16UL ) );  // View to a 16x16 submatrix of A

   // The function call operator provides access to all possible elements of the sparse submatrix,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse submatrix, the element is inserted into the
   // submatrix.
   sm(2,4) = 2.0;

   // An alternative for inserting elements into the submatrix is the \c insert() function. However,
   // it inserts the element only in case the element is not already contained in the submatrix.
   sm.insert( 2UL, 6UL, 3.7 );

   // As well as in the case of sparse matrices, elements can also be inserted via the \c append()
   // function. Also in case of submatrices, \c append() requires that the appended element's
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
// \c row() and \c column() functions, the current total capacity via the \c capacity() function,
// and the number of non-zero elements via the \c nonZeros() function. However, since submatrices
// are views to a specific submatrix of a matrix, several operations are not possible on views,
// such as resizing and swapping:

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>            SubmatrixType;

   MatrixType A;
   // ... Resizing and initialization

   // Creating a view to the a 8x12 submatrix of matrix A
   SubmatrixType sm = sub( A, 0UL, 0UL, 8UL, 12UL );

   sm.rows();      // Returns the number of rows of the submatrix
   sm.columns();   // Returns the number of columns of the submatrix
   sm.capacity();  // Returns the capacity of the submatrix
   sm.nonZeros();  // Returns the number of non-zero elements contained in the submatrix

   sm.resize( 10UL, 8UL );  // Compilation error: Cannot resize a submatrix of a matrix

   SubmatrixType sm2 = sub( A, 8UL, 0UL, 12UL, 8UL );
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
   SubmatrixType sm( sub( S1, 0UL, 0UL, 8UL, 8UL ) );  // View on the 8x8 submatrix of matrix S1
                                                       // starting from row 0 and column 0

   sub( S1, 0UL, 8UL, 8UL, 8UL ) = S2;  // Sparse matrix initialization of the 8x8 submatrix
                                        // starting in row 0 and column 8
   sm = D1;                             // Dense matrix initialization of the second 8x8 submatrix

   S3 = sm + S2;                              // Sparse matrix/sparse matrix addition
   D2 = D1  - sub( S1, 8UL, 0UL, 8UL, 8UL );  // Dense matrix/sparse matrix subtraction
   S2 = sm * sub( S1, 8UL, 8UL, 8UL, 8UL );   // Sparse matrix/sparse matrix multiplication

   sub( S1, 8UL, 0UL, 8UL, 8UL ) *= 2.0;     // In-place scaling of a submatrix of S1
   b = sub( S1, 8UL, 8UL, 8UL, 8UL ) * 2.0;  // Scaling of the a submatrix of S1
   b = 2.0 * sm;                             // Scaling of the a submatrix of S1

   sub( S1, 0UL, 8UL, 8UL, 8UL ) += S2;  // Addition assignment
   sub( S1, 8UL, 0UL, 8UL, 8UL ) -= D1;  // Subtraction assignment
   sub( S1, 8UL, 8UL, 8UL, 8UL ) *= sm;  // Multiplication assignment

   a = sub( S1, 4UL, 4UL, 8UL, 8UL ) * b;  // Sparse matrix/dense vector multiplication
   \endcode
*/
template< typename MT                                 // Type of the sparse matrix
        , bool SO = IsColumnMajorMatrix<MT>::value >  // Storage order
class SparseSubmatrix : public SparseMatrix< SparseSubmatrix<MT,SO>, SO >
                      , private Expression
{
 private:
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
   typedef SparseSubmatrix<MT,SO>              This;           //!< Type of this SparseSubmatrix instance.
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
   template< typename IteratorType  // Type of the sparse matrix iterator
           , bool ConstFlag >       // Constness flag
   class SubmatrixElement
   {
    public:
      //**Type definitions*************************************************************************
      typedef typename SelectType< ConstFlag, const ElementType&, ElementType& >::Type  ReferenceType;
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
   template< typename IteratorType  // Type of the sparse matrix iterator
           , bool ConstFlag >       // Constness flag
   class SubmatrixIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::forward_iterator_tag                 IteratorCategory;  //!< The iterator category.
      typedef SubmatrixElement<IteratorType,ConstFlag>  ValueType;         //!< Type of the underlying elements.
      typedef ValueType                                 PointerType;       //!< Pointer return type.
      typedef ValueType                                 ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                                 DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the SubmatrixIterator class.
      //
      // \param pos Iterator to the current sparse element.
      // \param offset The offset within the according row/column of the sparse matrix.
      */
      inline SubmatrixIterator( IteratorType pos, size_t offset )
         : pos_   ( pos    )  // Iterator to the current sparse element
         , offset_( offset )  // The offset of the according row/column of the sparse matrix
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename IteratorType2, bool ConstFlag2 >
      inline SubmatrixIterator( const SubmatrixIterator<IteratorType2,ConstFlag2>& it )
         : pos_   ( it.pos_    )  // Iterator to the current sparse element.
         , offset_( it.offset_ )  // The offset of the according row/column of the sparse matrix
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
      template< typename IteratorType2, bool ConstFlag2 >
      inline bool operator==( const SubmatrixIterator<IteratorType2,ConstFlag2>& rhs ) const {
         return pos_ == rhs.pos_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side submatrix iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename IteratorType2, bool ConstFlag2 >
      inline bool operator!=( const SubmatrixIterator<IteratorType2,ConstFlag2>& rhs ) const {
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

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;     //!< Iterator to the current sparse element.
      size_t       offset_;  //!< The offset of the according row/column of the sparse matrix.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      /*! \cond BLAZE_INTERNAL */
      template< typename IteratorType2, bool ConstFlag2 > friend class SubmatrixIterator;
      template< typename MT2, bool SO2 > friend class SparseSubmatrix;
      /*! \endcond */
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubmatrixIterator<typename MT::ConstIterator,true>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, SubmatrixIterator<typename MT::Iterator,false> >::Type  Iterator;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseSubmatrix( MT& matrix, size_t row, size_t column, size_t m, size_t n );
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
                              inline size_t           rows() const;
                              inline size_t           columns() const;
                              inline size_t           capacity() const;
                              inline size_t           capacity( size_t i ) const;
                              inline size_t           nonZeros() const;
                              inline size_t           nonZeros( size_t i ) const;
                              inline void             reset();
                              inline void             reset( size_t i );
                                     Iterator         insert ( size_t i, size_t j, const ElementType& value );
                              inline void             erase  ( size_t i, size_t j );
                              inline Iterator         erase  ( size_t i, Iterator pos );
                              inline Iterator         erase  ( size_t i, Iterator first, Iterator last );
                              inline void             reserve( size_t nonzeros );
                                     void             reserve( size_t i, size_t nonzeros );
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
   template< typename MT2, bool SO2 > inline void assign   ( const DenseMatrix<MT2,SO2>&  rhs );
   template< typename MT2, bool SO2 > inline void assign   ( const SparseMatrix<MT2,SO2>& rhs );
   template< typename MT2, bool SO2 > inline void addAssign( const DenseMatrix<MT2,SO2>&  rhs );
   template< typename MT2, bool SO2 > inline void addAssign( const SparseMatrix<MT2,SO2>& rhs );
   template< typename MT2, bool SO2 > inline void subAssign( const DenseMatrix<MT2,SO2>&  rhs );
   template< typename MT2, bool SO2 > inline void subAssign( const SparseMatrix<MT2,SO2>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT&          matrix_;  //!< The sparse matrix containing the submatrix.
   const size_t row_;     //!< The first row of the submatrix.
   const size_t column_;  //!< The first column of the submatrix.
   const size_t m_;       //!< The number of rows of the submatrix.
   const size_t n_;       //!< The number of columns of the submatrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE  ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
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
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline SparseSubmatrix<MT,SO>::SparseSubmatrix( MT& matrix, size_t row, size_t column, size_t m, size_t n )
   : matrix_( matrix )  // The sparse matrix containing the submatrix
   , row_   ( row    )  // The first row of the submatrix
   , column_( column )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
{
   if( m == 0UL || row    + m > matrix.rows() ||
       n == 0UL || column + n > matrix.columns() )
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Reference
   SparseSubmatrix<MT,SO>::operator()( size_t i, size_t j )
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstReference
   SparseSubmatrix<MT,SO>::operator()( size_t i, size_t j ) const
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator SparseSubmatrix<MT,SO>::begin( size_t i )
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );
      return Iterator( matrix_.lowerBound( i + row_, column_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid sparse submatrix column access index" );
      return Iterator( matrix_.lowerBound( row_, i + column_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator SparseSubmatrix<MT,SO>::begin( size_t i ) const
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );
      return Iterator( matrix_.lowerBound( i + row_, column_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid sparse submatrix column access index" );
      return Iterator( matrix_.lowerBound( row_, i + column_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator SparseSubmatrix<MT,SO>::cbegin( size_t i ) const
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );
      return Iterator( matrix_.lowerBound( i + row_, column_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid sparse submatrix column access index" );
      return Iterator( matrix_.lowerBound( row_, i + column_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator SparseSubmatrix<MT,SO>::end( size_t i )
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );
      return Iterator( matrix_.lowerBound( i + row_, column_ + n_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid sparse submatrix column access index" );
      return Iterator( matrix_.lowerBound( row_ + m_, i + column_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator SparseSubmatrix<MT,SO>::end( size_t i ) const
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );
      return ConstIterator( matrix_.lowerBound( i + row_, column_ + n_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid sparse submatrix column access index" );
      return ConstIterator( matrix_.lowerBound( row_ + m_, i + column_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator SparseSubmatrix<MT,SO>::cend( size_t i ) const
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid sparse submatrix row access index" );
      return ConstIterator( matrix_.lowerBound( i + row_, column_ + n_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid sparse submatrix column access index" );
      return ConstIterator( matrix_.lowerBound( row_ + m_, i + column_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline SparseSubmatrix<MT,SO>& SparseSubmatrix<MT,SO>::operator=( const SparseSubmatrix& rhs )
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
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The sparse submatrix is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,SO>& SparseSubmatrix<MT,SO>::operator=( const Matrix<MT2,SO2>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename MT2::ResultType );

   if( rows() != (~rhs).rows() || columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,SO>& SparseSubmatrix<MT,SO>::operator+=( const Matrix<MT2,SO2>& rhs )
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,SO>& SparseSubmatrix<MT,SO>::operator-=( const Matrix<MT2,SO2>& rhs )
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline SparseSubmatrix<MT,SO>& SparseSubmatrix<MT,SO>::operator*=( const Matrix<MT2,SO2>& rhs )
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,SO> >::Type&
   SparseSubmatrix<MT,SO>::operator*=( Other rhs )
{
   const size_t iend( ( SO == rowMajor )?( rows() ):( columns() ) );

   for( size_t i=0UL; i<iend; ++i ) {
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubmatrix<MT,SO> >::Type&
   SparseSubmatrix<MT,SO>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   const size_t iend( ( SO == rowMajor )?( rows() ):( columns() ) );

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( size_t i=0UL; i<iend; ++i ) {
         const Iterator last( end(i) );
         for( Iterator element=begin(i); element!=last; ++element )
            element->value() *= tmp;
      }
   }
   else {
      for( size_t i=0UL; i<iend; ++i ) {
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
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,SO>::rows() const
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
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,SO>::columns() const
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
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,SO>::capacity() const
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
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,SO>::capacity( size_t i ) const
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
      return columns();
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
      return rows();
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the sparse submatrix
//
// \return The number of non-zero elements in the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,SO>::nonZeros() const
{
   const size_t iend( ( SO == rowMajor )?( rows() ):( columns() ) );
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<iend; ++i )
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
        , bool SO >    // Storage order
inline size_t SparseSubmatrix<MT,SO>::nonZeros( size_t i ) const
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
   }

   return end(i) - begin(i);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,SO>::reset()
{
   if( SO == rowMajor ) {
      for( size_t i=row_; i<row_+m_; ++i ) {
         matrix_.erase( i, matrix_.lowerBound( i, column_ ), matrix_.lowerBound( i, column_+n_ ) );
      }
   }
   else {
      for( size_t j=column_; j<column_+n_; ++j ) {
         matrix_.erase( j, matrix_.lowerBound( row_, j ), matrix_.lowerBound( row_+m_, j ) );
      }
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
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,SO>::reset( size_t i )
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
      const size_t index( row_ + i );
      matrix_.erase( index, matrix_.lowerBound( index, column_ ), matrix_.lowerBound( index, column_+n_ ) );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
      const size_t index( column_ + i );
      matrix_.erase( index, matrix_.lowerBound( row_, index ), matrix_.lowerBound( row_+m_, index ) );
   }
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
        , bool SO >    // Storage order
typename SparseSubmatrix<MT,SO>::Iterator
   SparseSubmatrix<MT,SO>::insert( size_t i, size_t j, const ElementType& value )
{
   const size_t offset( ( SO == rowMajor )?( column_ ):( row_ ) );
   return Iterator( matrix_.insert( row_+i, column_+j, value ), offset );
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
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,SO>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( row_ + i, column_ + j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse submatrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator
   SparseSubmatrix<MT,SO>::erase( size_t i, Iterator pos )
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
      return Iterator( matrix_.erase( row_+i, pos.pos_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
      return Iterator( matrix_.erase( column_+i, pos.pos_ ), row_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the sparse submatrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of element from the sparse submatrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator
   SparseSubmatrix<MT,SO>::erase( size_t i, Iterator first, Iterator last )
{
   if( SO == rowMajor ) {
      BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
      return Iterator( matrix_.erase( row_+i, first.pos_, last.pos_ ), column_ );
   }
   else {
      BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
      return Iterator( matrix_.erase( column_+i, first.pos_, last.pos_ ), row_ );
   }
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
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,SO>::reserve( size_t nonzeros )
{
   return;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of a specific row/column of the sparse submatrix.
//
// \param i The row/column index of the new element \f$[0..M-1]\f$ or \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified row.
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
        , bool SO >    // Storage order
void SparseSubmatrix<MT,SO>::reserve( size_t i, size_t nonzeros )
{
   return;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the sparse submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the sparse submatrix.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline SparseSubmatrix<MT,SO>& SparseSubmatrix<MT,SO>::scale( Other scalar )
{
   const size_t iend( ( SO == rowMajor )?( rows() ):( columns() ) );
   for( size_t i=0UL; i<iend; ++i ) {
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator
   SparseSubmatrix<MT,SO>::find( size_t i, size_t j )
{
   if( SO == rowMajor ) {
      const typename MT::Iterator pos( matrix_.find( row_ + i, column_ + j ) );

      if( pos != matrix_.end( row_ + i ) )
         return Iterator( pos, column_ );
      else
         return end( i );
   }
   else {
      const typename MT::Iterator pos( matrix_.find( row_ + i, column_ + j ) );

      if( pos != matrix_.end( column_ + j ) )
         return Iterator( pos, row_ );
      else
         return end( j );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator
   SparseSubmatrix<MT,SO>::find( size_t i, size_t j ) const
{
   if( SO == rowMajor ) {
      const typename MT::ConstIterator pos( matrix_.find( row_ + i, column_ + j ) );

      if( pos != matrix_.end( row_ + i ) )
         return ConstIterator( pos, column_ );
      else
         return end( i );
   }
   else {
      const typename MT::ConstIterator pos( matrix_.find( row_ + i, column_ + j ) );

      if( pos != matrix_.end( column_ + j ) )
         return ConstIterator( pos, row_ );
      else
         return end( j );
   }
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator
   SparseSubmatrix<MT,SO>::lowerBound( size_t i, size_t j )
{
   const size_t offset( ( SO == rowMajor )?( column_ ):( row_ ) );
   return Iterator( matrix_.lowerBound( row_ + i, column_ + j ), offset );
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator
   SparseSubmatrix<MT,SO>::lowerBound( size_t i, size_t j ) const
{
   const size_t offset( ( SO == rowMajor )?( column_ ):( row_ ) );
   return ConstIterator( matrix_.lowerBound( row_ + i, column_ + j ), offset );
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::Iterator
   SparseSubmatrix<MT,SO>::upperBound( size_t i, size_t j )
{
   const size_t offset( ( SO == rowMajor )?( column_ ):( row_ ) );
   return Iterator( matrix_.upperBound( row_ + i, column_ + j ), offset );
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
        , bool SO >    // Storage order
inline typename SparseSubmatrix<MT,SO>::ConstIterator
   SparseSubmatrix<MT,SO>::upperBound( size_t i, size_t j ) const
{
   const size_t offset( ( SO == rowMajor )?( column_ ):( row_ ) );
   return ConstIterator( matrix_.upperBound( row_ + i, column_ + j ), offset );
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
   SubmatrixType B = sub( A, 10, 10, 4, 3 );

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
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,SO>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault( value ) )
      matrix_.insert( row_ + i, column_ + j, value );
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
        , bool SO >    // Storage order
inline void SparseSubmatrix<MT,SO>::finalize( size_t i )
{
   return;
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,SO>::canAlias( const Other* alias ) const
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubmatrix<MT,SO>::isAliased( const Other* alias ) const
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,SO>::assign( const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         for( size_t j=0UL; j<columns(); ++j ) {
            append( i, j, (~rhs)(i,j), true );
         }
      }
   }
   else {
      for( size_t j=0UL; j<columns(); ++j ) {
         for( size_t i=0UL; i<rows(); ++i ) {
            append( i, j, (~rhs)(i,j), true );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT   // Type of the sparse matrix
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,SO>::assign( const SparseMatrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (~rhs).columns(), "Invalid number of columns" );

   if( SO2 == rowMajor ) {
      for( size_t i=0UL; i<(~rhs).rows(); ++i ) {
         for( typename MT2::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element ) {
            append( i, element->index(), element->value(), true );
         }
      }
   }
   else {
      for( size_t j=0UL; j<(~rhs).columns(); ++j ) {
         for( typename MT2::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element ) {
            append( element->index(), j, element->value(), true );
         }
      }
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,SO>::addAssign( const DenseMatrix<MT2,SO2>& rhs )
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand side sparse matrix
inline void SparseSubmatrix<MT,SO>::addAssign( const SparseMatrix<MT2,SO2>& rhs )
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline void SparseSubmatrix<MT,SO>::subAssign( const DenseMatrix<MT2,SO2>& rhs )
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
        , bool SO >     // Storage order
template< typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand sparse matrix
inline void SparseSubmatrix<MT,SO>::subAssign( const SparseMatrix<MT2,SO2>& rhs )
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
//  SPARSESUBMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseSubmatrix operators */
//@{
template< typename MT, bool SO >
inline void reset( SparseSubmatrix<MT,SO>& sm );

template< typename MT, bool SO >
inline void clear( SparseSubmatrix<MT,SO>& sm );

template< typename MT, bool SO >
inline bool isDefault( const SparseSubmatrix<MT,SO>& sm );
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
        , bool SO >    // Storage order
inline void reset( SparseSubmatrix<MT,SO>& sm )
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
        , bool SO >    // Storage order
inline void clear( SparseSubmatrix<MT,SO>& sm )
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
   if( isDefault( sub( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline bool isDefault( const SparseSubmatrix<MT,SO>& sm )
{
   using blaze::isDefault;

   typedef typename SparseSubmatrix<MT,SO>::ConstIterator  ConstIterator;

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
//  GLOBAL FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given sparse matrix.
// \ingroup views
//
// \param sm The sparse matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specific submatrix of the sparse matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given sparse
// matrix. The following example demonstrates the creation of a submatrix of size 4 by 4 starting
// from position (3,2):

   \code
   using blaze::rowMatrix;

   typedef blaze::CompressedMatrix<double,rowMatrix>  Matrix;

   Matrix A;
   // ... Resizing and initialization
   blaze::SparseSubmatrix<Matrix> = sub( A, 3UL, 2UL, 4UL, 4UL );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename DisableIf< Or< IsComputation<MT>, IsTransExpr<MT> >, SparseSubmatrix<MT> >::Type
   sub( SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return SparseSubmatrix<MT>( ~sm, row, column, m, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given sparse matrix.
// \ingroup views
//
// \param sm The sparse matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specific submatrix of the sparse matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given sparse
// matrix. The following example demonstrates the creation of a submatrix of size 4 by 4 starting
// from position (3,2):

   \code
   using blaze::rowMatrix;

   typedef blaze::CompressedMatrix<double,rowMatrix>  Matrix;

   Matrix A;
   // ... Resizing and initialization
   blaze::SparseSubmatrix<Matrix> = sub( A, 3UL, 2UL, 4UL, 4UL );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename DisableIf< Or< IsComputation<MT>, IsTransExpr<MT> >, SparseSubmatrix<const MT> >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return SparseSubmatrix<const MT>( ~sm, row, column, m, n );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix addition.
// \ingroup views
//
// \param sm The constant matrix/matrix addition.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the addition.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix addition.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatAddExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return sub( (~sm).leftOperand(), row, column, m, n ) + sub( (~sm).rightOperand(), row, column, m, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix subtraction.
// \ingroup views
//
// \param sm The constant matrix/matrix subtraction.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the subtraction.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix subtraction.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatSubExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return sub( (~sm).leftOperand(), row, column, m, n ) - sub( (~sm).rightOperand(), row, column, m, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix multiplication.
// \ingroup views
//
// \param sm The constant matrix/matrix multiplication.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the multiplication.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix multiplication.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatMultExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   typename MT::LeftOperand  left ( (~sm).leftOperand()  );
   typename MT::RightOperand right( (~sm).rightOperand() );

   return sub( left, row, 0UL, m, left.columns() ) * sub( right, 0UL, column, right.rows(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given outer product.
// \ingroup views
//
// \param sm The constant outer product.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the outer product.
//
// This function returns an expression representing the specified submatrix of the given
// outer product.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsVecTVecMultExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return sub( (~sm).leftOperand(), row, m ) * sub( (~sm).rightOperand(), column, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/scalar multiplication.
// \ingroup views
//
// \param sm The constant matrix/scalar multiplication.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the multiplication.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/scalar multiplication.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarMultExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return sub( (~sm).leftOperand(), row, column, m, n ) * (~sm).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/scalar division.
// \ingroup views
//
// \param sm The constant matrix/scalar division.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the division.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/scalar division.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarDivExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return sub( (~sm).leftOperand(), row, column, m, n ) / (~sm).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix abs operation.
// \ingroup views
//
// \param sv The constant matrix abs operation.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the abs operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// abs operation.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatAbsExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return abs( sub( (~sm).operand(), row, column, m, n ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix evaluation operation.
// \ingroup views
//
// \param sv The constant matrix evaluation operation.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the evaluation operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// evaluation operation.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatEvalExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return eval( sub( (~sm).operand(), row, column, m, n ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix transpose operation.
// \ingroup views
//
// \param sv The constant matrix transpose operation.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the transpose operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// transpose operation.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatTransExpr<MT>, typename SubmatrixExprTrait<MT>::Type >::Type
   sub( const SparseMatrix<MT,SO>& sm, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return trans( sub( (~sm).operand(), column, row, n, m ) );
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
template< typename MT, bool SO >
struct SubmatrixTrait< SparseSubmatrix<MT,SO> >
{
   typedef typename SubvectorTrait< typename SparseSubmatrix<MT,SO>::ResultType >::Type  Type;
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
template< typename MT, bool SO >
struct RowTrait< SparseSubmatrix<MT,SO> >
{
   typedef typename RowTrait< typename SparseSubmatrix<MT,SO>::ResultType >::Type  Type;
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
template< typename MT, bool SO >
struct ColumnTrait< SparseSubmatrix<MT,SO> >
{
   typedef typename ColumnTrait< typename SparseSubmatrix<MT,SO>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
