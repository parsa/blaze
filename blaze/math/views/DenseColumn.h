//=================================================================================================
/*!
//  \file blaze/math/views/DenseColumn.h
//  \brief Header file for the DenseColumn class template
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

#ifndef _BLAZE_MATH_VIEWS_DENSECOLUMN_H_
#define _BLAZE_MATH_VIEWS_DENSECOLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <stdexcept>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/Forward.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
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
/*!\defgroup dense_column DenseColumn
// \ingroup views
*/
/*!\brief Reference to a specific column of a dense matrix.
// \ingroup dense_column
//
// The DenseRow template represents a reference to a specific column of a dense matrix primitive.
// The type of the dense matrix is specified via the first template parameter:

   \code
   template< typename MT, bool SO >
   class DenseColumn;
   \endcode

//  - MT: specifies the type of the dense matrix primitive. DenseColumn can be used with any
//        dense matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the dense matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//
// A reference to a dense column can be conventiently created via the column() function. The
// column can be either used as an alias to grant write access to a specific column of a matrix
// primitive on the left-hand side of an assignment or to grant read-access to a specific column
// of a matrix primitive or expression on the right-hand side of an assignment:

   \code
   blaze::DynamicVector<double,columnVector> x;
   blaze::DynamicMatrix<double,columnMajor> A, B;
   // ... Resizing and initialization

   // Setting the 2nd column of matrix A to x
   column( A, 2UL ) = x;

   // Setting x to the 3rd column of the result of the matrix multiplication
   x = column( A * B, 3UL );
   \endcode

// A dense column can be used like any other column vector. The elements of the dense column can
// be directly accessed with the subscript operator. The numbering of the column elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of rows of the referenced matrix. The following example gives an
// impression of the use of DenseColumn. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse vectors with
// fitting element types:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::DynamicMatrix;

   DynamicVector<double,columnVector> a( 2UL, 2.0 ), b;
   CompressedVector<double,columnVector> c( 2UL );
   c[1] = 3.0;

   typedef DynamicMatrix<double,columnMajor>  MatrixType;
   MatrixType A( 2UL, 4UL );  // Non-initialized 2x4 matrix

   typedef DenseColumn<DenseMatrix>  RowType;
   RowType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   col0[0] = 0UL;           // Manual initialization of the 0th column of A
   col0[1] = 0UL;
   column( A, 1UL ) = 1.0;  // Homogeneous initialization of the 1st column of A
   column( A, 2UL ) = a;    // Dense vector initialization of the 2nd column of A
   column( A, 3UL ) = c;    // Sparse vector initialization of the 3rd column of A

   b = col0 + a;                 // Dense vector/dense vector addition
   b = c + column( A, 1UL );     // Sparse vector/dense vector addition
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

// It is possible to create a column view on both row-major and column-major matrices. However,
// please note that creating a column view on a matrix stored in row-major fashion can result
// in a considerable performance decrease in comparison to a column view on a column-major matrix
// due to the non-contiguous storage of the matrix elements. Therefore care has to be taken in
// the choice of the most suitable storage order:

   \code
   // Setup of two row-major matrices
   DynamicMatrix<double,rowMajor> A( 128UL, 128UL );
   DynamicMatrix<double,rowMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th column of the multiplication between A and B ...
   DynamicVector<double,columnVector> x = column( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // A with the 15th column of the row-major matrix B.
   DynamicVector<double,rowVector> x = A * column( B, 15UL );
   \endcode

// Although Blaze performs the resulting matrix/vector multiplication as efficiently as possible
// using a column-major storage order for matrix B would result in a more efficient evaluation.
*/
template< typename MT                                 // Type of the dense matrix
        , bool SO = IsColumnMajorMatrix<MT>::value >  // Storage order
class DenseColumn : public DenseVector< DenseColumn<MT,SO>, false >
                  , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef IntrinsicTrait<typename MT::ElementType>  IT;  //!< Intrinsic trait for the column element type.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense column will
       return references and iterators to const. Otherwise \a useConst will be set to 0
       and the dense column will offer write access to the dense matrix elements both via
       the subscript operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseColumn<MT,SO>                  This;           //!< Type of this DenseColumn instance.
   typedef typename ColumnTrait<MT>::Type      ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the column elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the column elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseColumn&                  CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant column value.
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, typename MT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = MT::vectorizable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseColumn( MT& matrix, size_t index );
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
                           inline DenseColumn& operator= ( const ElementType& rhs );
                           inline DenseColumn& operator= ( const DenseColumn& rhs );
   template< typename VT > inline DenseColumn& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline DenseColumn& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline DenseColumn& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline DenseColumn& operator*=( const Vector<VT,false>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseColumn >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseColumn >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t       size() const;
                              inline size_t       capacity() const;
                              inline size_t       nonZeros() const;
                              inline void         reset();
   template< typename Other > inline DenseColumn& scale( const Other& scalar );
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
   template< typename Other > inline bool          canAlias ( const Other* alias ) const;
   template< typename Other > inline bool          isAliased( const Other* alias ) const;
                              inline IntrinsicType get      ( size_t index ) const;

   template< typename VT >
   inline typename DisableIf< VectorizedAssign<VT> >::Type
      assign( const DenseVector<VT,false>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedAssign<VT> >::Type
      assign( const DenseVector<VT,false>& rhs );

   template< typename VT > inline void assign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedAddAssign<VT> >::Type
      addAssign( const DenseVector<VT,false>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedAddAssign<VT> >::Type
      addAssign( const DenseVector<VT,false>& rhs );

   template< typename VT > inline void addAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedSubAssign<VT> >::Type
      subAssign( const DenseVector<VT,false>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedSubAssign<VT> >::Type
      subAssign( const DenseVector<VT,false>& rhs );

   template< typename VT > inline void subAssign( const SparseVector<VT,false>& rhs );

   template< typename VT >
   inline typename DisableIf< VectorizedMultAssign<VT> >::Type
      multAssign( const DenseVector<VT,false>& rhs );

   template< typename VT >
   inline typename EnableIf< VectorizedMultAssign<VT> >::Type
      multAssign( const DenseVector<VT,false>& rhs );

   template< typename VT > inline void multAssign( const SparseVector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT&          matrix_;  //!< The dense matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE     ( MT );
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
/*!\brief The constructor for DenseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseColumn<MT,SO>::DenseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the column
   , col_   ( index  )  // The index of the column in the matrix
{
   if( matrix_.columns() <= index )
      throw std::invalid_argument( "Invalid column access index" );
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
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::Reference DenseColumn<MT,SO>::operator[]( size_t index )
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
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::ConstReference DenseColumn<MT,SO>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,col_);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::Iterator DenseColumn<MT,SO>::begin()
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
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::ConstIterator DenseColumn<MT,SO>::begin() const
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
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::ConstIterator DenseColumn<MT,SO>::cbegin() const
{
   return matrix_.begin( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::Iterator DenseColumn<MT,SO>::end()
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
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::ConstIterator DenseColumn<MT,SO>::end() const
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
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::ConstIterator DenseColumn<MT,SO>::cend() const
{
   return matrix_.end( col_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Homogenous assignment to all column elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::operator=( const ElementType& rhs )
{
   const size_t rows( size() );

   for( size_t i=0UL; i<rows; ++i )
      matrix_(i,col_) = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DenseColumn.
//
// \param rhs Dense column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::operator=( const DenseColumn& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Column sizes do not match" );

   const size_t rows( size() );

   for( size_t i=0UL; i<rows; ++i )
      matrix_(i,col_) = rhs[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::operator=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      const typename VT::ResultType tmp( ~rhs );
      assign( *this, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::operator+=( const Vector<VT,false>& rhs )
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      const typename VT::ResultType tmp( ~rhs );
      addAssign( *this, tmp );
   }
   else {
      addAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::operator-=( const Vector<VT,false>& rhs )
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      const typename VT::ResultType tmp( ~rhs );
      subAssign( *this, tmp );
   }
   else {
      subAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the dense matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::multAssign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) || IsSparseVector<VT>::value ) {
      const typename VT::ResultType tmp( ~rhs );
      multAssign( *this, tmp );
   }
   else {
      multAssign( *this, ~rhs );
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
inline typename EnableIf< IsNumeric<Other>, DenseColumn<MT,SO> >::Type&
   DenseColumn<MT,SO>::operator*=( Other rhs )
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
inline typename EnableIf< IsNumeric<Other>, DenseColumn<MT,SO> >::Type&
   DenseColumn<MT,SO>::operator/=( Other rhs )
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
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline size_t DenseColumn<MT,SO>::size() const
{
   return matrix_.rows();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the dense column.
//
// \return The capacity of the dense column.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline size_t DenseColumn<MT,SO>::capacity() const
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
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline size_t DenseColumn<MT,SO>::nonZeros() const
{
   return matrix_.nonZeros( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline void DenseColumn<MT,SO>::reset()
{
   matrix_.reset( col_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the dense column.
*/
template< typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline DenseColumn<MT,SO>& DenseColumn<MT,SO>::scale( const Other& scalar )
{
   for( size_t j=0UL; j<size(); ++j ) {
      matrix_(j,col_) *= scalar;
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
/*!\brief Returns whether the dense column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In
// contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DenseColumn<MT,SO>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In
// contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DenseColumn<MT,SO>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the intrinsic elements of the dense column.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed values.
//
// This function offers a direct access to the intrinsic elements of the column. It must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DenseColumn<MT,SO>::IntrinsicType DenseColumn<MT,SO>::get( size_t index ) const
{
   return matrix_.get( index, col_ );
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
inline typename DisableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseColumn<MT,SO>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) = (~rhs)[i    ];
      matrix_(i+1UL,col_) = (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) = (~rhs)[iend];
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
inline typename EnableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedAssign<VT> >::Type
   DenseColumn<MT,SO>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t rows( size() );

   if( rows > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(~rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<rows; i+=IT::size ) {
         stream( &matrix_(i,col_), (~rhs).get(i) );
      }
   }
   else
   {
      BLAZE_INTERNAL_ASSERT( ( rows - ( rows % (IT::size*4UL) ) ) == ( rows & size_t(-IT::size*4) ), "Invalid end calculation" );
      const size_t iend( rows & size_t(-IT::size*4) );

      for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
         store( &matrix_(i             ,col_), (~rhs).get(i             ) );
         store( &matrix_(i+IT::size    ,col_), (~rhs).get(i+IT::size    ) );
         store( &matrix_(i+IT::size*2UL,col_), (~rhs).get(i+IT::size*2UL) );
         store( &matrix_(i+IT::size*3UL,col_), (~rhs).get(i+IT::size*3UL) );
      }
      for( size_t i=iend; i<rows; i+=IT::size ) {
         store( &matrix_(i,col_), (~rhs).get(i) );
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
inline void DenseColumn<MT,SO>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) = element->value();
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
inline typename DisableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseColumn<MT,SO>::addAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) += (~rhs)[i    ];
      matrix_(i+1UL,col_) += (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) += (~rhs)[iend];
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
inline typename EnableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedAddAssign<VT> >::Type
   DenseColumn<MT,SO>::addAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t rows( size() );

   BLAZE_INTERNAL_ASSERT( ( rows - ( rows % (IT::size*4UL) ) ) == ( rows & size_t(-IT::size*4) ), "Invalid end calculation" );
   const size_t iend( rows & size_t(-IT::size*4) );

   for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
      store( &matrix_(i             ,col_), load( &matrix_(i             ,col_) ) + (~rhs).get(i             ) );
      store( &matrix_(i+IT::size    ,col_), load( &matrix_(i+IT::size    ,col_) ) + (~rhs).get(i+IT::size    ) );
      store( &matrix_(i+IT::size*2UL,col_), load( &matrix_(i+IT::size*2UL,col_) ) + (~rhs).get(i+IT::size*2UL) );
      store( &matrix_(i+IT::size*3UL,col_), load( &matrix_(i+IT::size*3UL,col_) ) + (~rhs).get(i+IT::size*3UL) );
   }
   for( size_t i=iend; i<rows; i+=IT::size ) {
      store( &matrix_(i,col_), load( &matrix_(i,col_) ) + (~rhs).get(i) );
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
inline void DenseColumn<MT,SO>::addAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) += element->value();
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
inline typename DisableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseColumn<MT,SO>::subAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) -= (~rhs)[i    ];
      matrix_(i+1UL,col_) -= (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) -= (~rhs)[iend];
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
inline typename EnableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedSubAssign<VT> >::Type
   DenseColumn<MT,SO>::subAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t rows( size() );

   BLAZE_INTERNAL_ASSERT( ( rows - ( rows % (IT::size*4UL) ) ) == ( rows & size_t(-IT::size*4) ), "Invalid end calculation" );
   const size_t iend( rows & size_t(-IT::size*4) );

   for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
      store( &matrix_(i             ,col_), load( &matrix_(i             ,col_) ) - (~rhs).get(i             ) );
      store( &matrix_(i+IT::size    ,col_), load( &matrix_(i+IT::size    ,col_) ) - (~rhs).get(i+IT::size    ) );
      store( &matrix_(i+IT::size*2UL,col_), load( &matrix_(i+IT::size*2UL,col_) ) - (~rhs).get(i+IT::size*2UL) );
      store( &matrix_(i+IT::size*3UL,col_), load( &matrix_(i+IT::size*3UL,col_) ) - (~rhs).get(i+IT::size*3UL) );
   }
   for( size_t i=iend; i<rows; i+=IT::size ) {
      store( &matrix_(i,col_), load( &matrix_(i,col_) ) - (~rhs).get(i) );
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
inline void DenseColumn<MT,SO>::subAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) -= element->value();
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
inline typename DisableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseColumn<MT,SO>::multAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) *= (~rhs)[i    ];
      matrix_(i+1UL,col_) *= (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) *= (~rhs)[iend];
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
inline typename EnableIf< typename DenseColumn<MT,SO>::BLAZE_TEMPLATE VectorizedMultAssign<VT> >::Type
   DenseColumn<MT,SO>::multAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t rows( size() );

   BLAZE_INTERNAL_ASSERT( ( rows - ( rows % (IT::size*4UL) ) ) == ( rows & size_t(-IT::size*4) ), "Invalid end calculation" );
   const size_t iend( rows & size_t(-IT::size*4) );

   for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
      store( &matrix_(i             ,col_), load( &matrix_(i             ,col_) ) * (~rhs).get(i             ) );
      store( &matrix_(i+IT::size    ,col_), load( &matrix_(i+IT::size    ,col_) ) * (~rhs).get(i+IT::size    ) );
      store( &matrix_(i+IT::size*2UL,col_), load( &matrix_(i+IT::size*2UL,col_) ) * (~rhs).get(i+IT::size*2UL) );
      store( &matrix_(i+IT::size*3UL,col_), load( &matrix_(i+IT::size*3UL,col_) ) * (~rhs).get(i+IT::size*3UL) );
   }
   for( size_t i=iend; i<rows; i+=IT::size ) {
      store( &matrix_(i,col_), load( &matrix_(i,col_) ) * (~rhs).get(i) );
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
inline void DenseColumn<MT,SO>::multAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( *this );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) = tmp[element->index()] * element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseColumn for row-major matrices.
// \ingroup views
//
// This specialization of DenseColumn adapts the class template to the requirements of
// row-major matrices.
*/
template< typename MT >  // Type of the dense matrix
class DenseColumn<MT,false> : public DenseVector< DenseColumn<MT,false>, false >
                            , private Expression
{
 private:
   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the dense column will
       return references and iterators to const. Otherwise \a useConst will be set to 0
       and the dense column will offer write access to the dense matrix elements both via
       the subscript operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseColumn<MT,false>               This;            //!< Type of this DenseColumn instance.
   typedef typename ColumnTrait<MT>::Type      ResultType;      //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;     //!< Type of the column elements.
   typedef typename MT::ReturnType             ReturnType;      //!< Return type for expression template evaluations
   typedef const ResultType                    CompositeType;   //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant column value.
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;
   //**********************************************************************************************

   //**ColumnIterator class definition*************************************************************
   /*!\brief Iterator over the elements of the dense column.
   */
   template< typename MatrixType >  // Type of the dense matrix
   class ColumnIterator
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

      typedef std::forward_iterator_tag   IteratorCategory;  //!< The iterator category.
      typedef RemoveReference<Reference>  ValueType;         //!< Type of the underlying elements.
      typedef ValueType*                  PointerType;       //!< Pointer return type.
      typedef Reference                   ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                   DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( matrix )  // The dense matrix containing the column.
         , row_   ( row    )  // The current row index.
         , column_( column )  // The current column index.
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ColumnIterator instances.
      //
      // \param it The column iterator to be copied.
      */
      template< typename MatrixType2 >
      inline ColumnIterator( const ColumnIterator<MatrixType2>& it )
         : matrix_( it.matrix_ )  // The dense matrix containing the column.
         , row_   ( it.row_    )  // The current row index.
         , column_( it.column_ )  // The current column index.
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ColumnIterator& operator++() {
         ++row_;
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
      /*!\brief Direct access to the dense vector element at the current iterator position.
      //
      // \return The current value of the dense element.
      */
      inline ReferenceType operator*() const {
         return matrix_(row_,column_);
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the dense vector element at the current iterator position.
      //
      // \return Reference to the dense vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return &matrix_(row_,column_);
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2 >
      inline bool operator==( const ColumnIterator<MatrixType2>& rhs ) const {
         return ( &matrix_ == &rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ == rhs.column_ );
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2 >
      inline bool operator!=( const ColumnIterator<MatrixType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two expression iterators.
      //
      // \param rhs The right-hand side expression iterator.
      // \return The number of elements between the two expression iterators.
      */
      inline DifferenceType operator-( const ColumnIterator& rhs ) const {
         return row_ - rhs.row_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType&  matrix_;  //!< The dense matrix containing the column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      /*! \cond BLAZE_INTERNAL */
      template< typename MatrixType2 > friend class ColumnIterator;
      /*! \endcond */
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef ColumnIterator<const MT>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, ColumnIterator<MT> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseColumn( MT& matrix, size_t index );
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
                           inline DenseColumn& operator= ( const ElementType& rhs );
                           inline DenseColumn& operator= ( const DenseColumn& rhs );
   template< typename VT > inline DenseColumn& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline DenseColumn& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline DenseColumn& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline DenseColumn& operator*=( const Vector<VT,false>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseColumn >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseColumn >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t       size() const;
                              inline size_t       capacity() const;
                              inline size_t       nonZeros() const;
                              inline void         reset();
   template< typename Other > inline DenseColumn& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias  ( const Other* alias ) const;
   template< typename Other > inline bool isAliased ( const Other* alias ) const;
   template< typename VT >    inline void assign    ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign    ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void addAssign ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void subAssign ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void subAssign ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void multAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void multAssign( const SparseVector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT&          matrix_;  //!< The dense matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE  ( MT );
   /*! \endcond */
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
/*!\brief The constructor for DenseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT >  // Type of the dense matrix
inline DenseColumn<MT,false>::DenseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The dense matrix containing the column
   , col_   ( index  )  // The index of the column in the matrix
{
   if( matrix_.columns() <= index )
      throw std::invalid_argument( "Invalid column access index" );
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
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::Reference DenseColumn<MT,false>::operator[]( size_t index )
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
*/
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::ConstReference DenseColumn<MT,false>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,col_);
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
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::Iterator DenseColumn<MT,false>::begin()
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
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::ConstIterator DenseColumn<MT,false>::begin() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::ConstIterator DenseColumn<MT,false>::cbegin() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::Iterator DenseColumn<MT,false>::end()
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
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::ConstIterator DenseColumn<MT,false>::end() const
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
template< typename MT >  // Type of the dense matrix
inline typename DenseColumn<MT,false>::ConstIterator DenseColumn<MT,false>::cend() const
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
/*!\brief Homogenous assignment to all column elements.
//
// \param rhs Scalar value to be assigned to all column elements.
// \return Reference to the assigned column.
*/
template< typename MT >  // Type of the dense matrix
inline DenseColumn<MT,false>& DenseColumn<MT,false>::operator=( const ElementType& rhs )
{
   const size_t rows( size() );

   for( size_t i=0UL; i<rows; ++i )
      matrix_(i,col_) = rhs;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseColumn.
//
// \param rhs Dense column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
inline DenseColumn<MT,false>& DenseColumn<MT,false>::operator=( const DenseColumn& rhs )
{
   if( &rhs == this ) return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Column sizes do not match" );

   const size_t rows( size() );

   for( size_t i=0UL; i<rows; ++i )
      matrix_(i,col_) = rhs[i];

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,false>& DenseColumn<MT,false>::operator=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      const ResultType tmp( ~rhs );
      assign( *this, tmp );
   }
   else {
      if( IsSparseVector<VT>::value )
         reset();
      assign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,false>& DenseColumn<MT,false>::operator+=( const Vector<VT,false>& rhs )
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      const typename VT::ResultType tmp( ~rhs );
      addAssign( *this, tmp );
   }
   else {
      addAssign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,false>& DenseColumn<MT,false>::operator-=( const Vector<VT,false>& rhs )
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      const typename VT::ResultType tmp( ~rhs );
      subAssign( *this, tmp );
   }
   else {
      subAssign( *this, ~rhs );
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
// \param rhs The right-hand side vector to be multiplied with the dense column.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the dense matrix
template< typename VT >  // Type of the right-hand side vector
inline DenseColumn<MT,false>& DenseColumn<MT,false>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::multAssign;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( typename VT::ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( this ) || IsSparseVector<VT>::value ) {
      const typename VT::ResultType tmp( ~rhs );
      multAssign( *this, tmp );
   }
   else {
      multAssign( *this, ~rhs );
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
inline typename EnableIf< IsNumeric<Other>, DenseColumn<MT,false> >::Type&
   DenseColumn<MT,false>::operator*=( Other rhs )
{
   return operator=( (*this) * rhs );
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
inline typename EnableIf< IsNumeric<Other>, DenseColumn<MT,false> >::Type&
   DenseColumn<MT,false>::operator/=( Other rhs )
{
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
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseColumn<MT,false>::size() const
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense column.
//
// \return The capacity of the dense column.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseColumn<MT,false>::capacity() const
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
// of columns of the matrix containing the column.
*/
template< typename MT >  // Type of the dense matrix
inline size_t DenseColumn<MT,false>::nonZeros() const
{
   const size_t rows( size() );
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows; ++i )
      if( !isDefault( matrix_(i,col_) ) )
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
inline void DenseColumn<MT,false>::reset()
{
   using blaze::reset;
   const size_t rows( size() );
   for( size_t i=0UL; i<rows; ++i )
      reset( matrix_(i,col_) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the dense column.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the scalar value
inline DenseColumn<MT,false>& DenseColumn<MT,false>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<size(); ++i ) {
      matrix_(i,col_) *= scalar;
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
/*!\brief Returns whether the dense column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address can alias with the dense column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseColumn<MT,false>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense column, \a false if not.
//
// This function returns whether the given address is aliased with the dense column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT >     // Type of the dense matrix
template< typename Other >  // Data type of the foreign expression
inline bool DenseColumn<MT,false>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
inline void DenseColumn<MT,false>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) = (~rhs)[i    ];
      matrix_(i+1UL,col_) = (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) = (~rhs)[iend];
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
inline void DenseColumn<MT,false>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) = element->value();
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
inline void DenseColumn<MT,false>::addAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) += (~rhs)[i    ];
      matrix_(i+1UL,col_) += (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) += (~rhs)[iend];
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
inline void DenseColumn<MT,false>::addAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) += element->value();
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
inline void DenseColumn<MT,false>::subAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) -= (~rhs)[i    ];
      matrix_(i+1UL,col_) -= (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) -= (~rhs)[iend];
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
inline void DenseColumn<MT,false>::subAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) -= element->value();
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
inline void DenseColumn<MT,false>::multAssign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( (~rhs).size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      matrix_(i    ,col_) *= (~rhs)[i    ];
      matrix_(i+1UL,col_) *= (~rhs)[i+1UL];
   }
   if( iend < (~rhs).size() )
      matrix_(iend,col_) *= (~rhs)[iend];
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
inline void DenseColumn<MT,false>::multAssign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( *this );

   reset();

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      matrix_(element->index(),col_) = tmp[element->index()] * element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  DENSECOLUMN OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseColumn operators */
//@{
template< typename MT, bool SO >
inline void reset( DenseColumn<MT,SO>& column );

template< typename MT, bool SO >
inline void clear( DenseColumn<MT,SO>& column );

template< typename MT, bool SO >
inline bool isnan( const DenseColumn<MT,SO>& column );

template< typename MT, bool SO >
inline bool isDefault( const DenseColumn<MT,SO>& column );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dense column.
// \ingroup dense_column
//
// \param column The dense column to be resetted.
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline void reset( DenseColumn<MT,SO>& column )
{
   column.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dense column.
// \ingroup dense_column
//
// \param column The dense column to be cleared.
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline void clear( DenseColumn<MT,SO>& column )
{
   column.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense column for not-a-number elements.
// \ingroup dense_column
//
// \param column The dense column to be checked for not-a-number elements.
// \return \a true if at least one element of the column is not-a-number, \a false otherwise.
//
// This function checks the dense column for not-a-number (NaN) elements. If at least one element
// of the column is not-a-number, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isnan( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline bool isnan( const DenseColumn<MT,SO>& column )
{
   for( size_t i=0UL; i<column.size(); ++i ) {
      if( isnan( column[i] ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense column is in default state.
// \ingroup dense_column
//
// \param column The dense column to be tested for its default state.
// \return \a true in case the given dense column is component-wise zero, \a false otherwise.
//
// This function checks whether the dense column is in default state. For instance, in case the
// column is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all column elements are 0 and \a false in case any column element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline bool isDefault( const DenseColumn<MT,SO>& column )
{
   for( size_t i=0UL; i<column.size(); ++i )
      if( !isDefault( column[i] ) ) return false;
   return true;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given dense matrix.
// \ingroup views
//
// \param dm The dense matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given dense matrix.

   \code
   typedef blaze::DynamicMatrix<double,columnMajor>  Matrix;

   Matrix A;
   // ... Resizing and initialization
   blaze::DenseColumn<Matrix> = column( A, 3UL );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseColumn<MT> column( DenseMatrix<MT,SO>& dm, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return DenseColumn<MT>( ~dm, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant dense matrix.
// \ingroup views
//
// \param dm The constant dense matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given dense matrix.

   \code
   typedef blaze::DynamicMatrix<double,columnMajor>  Matrix;

   const Matrix A;
   // ... Resizing and initialization
   blaze::DenseColumn<const Matrix> = column( A, 3UL );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline DenseColumn<const MT> column( const DenseMatrix<MT,SO>& dm, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return DenseColumn<const MT>( ~dm, index );
}
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2, size_t N >
struct AddTrait< DenseColumn<T1,SO>, StaticVector<T2,N,false> >
{
   typedef typename AddTrait< typename DenseColumn<T1,SO>::ResultType,
                              StaticVector<T2,N,false> >::Type  Type;
};

template< typename T1, size_t N, typename T2, bool SO >
struct AddTrait< StaticVector<T1,N,false>, DenseColumn<T2,SO> >
{
   typedef typename AddTrait< StaticVector<T1,N,false>,
                              typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< DenseColumn<T1,SO>, DynamicVector<T2,false> >
{
   typedef typename AddTrait< typename DenseColumn<T1,SO>::ResultType,
                              DynamicVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct AddTrait< DynamicVector<T1,false>, DenseColumn<T2,SO> >
{
   typedef typename AddTrait< DynamicVector<T1,false>,
                              typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< DenseColumn<T1,SO>, CompressedVector<T2,false> >
{
   typedef typename AddTrait< typename DenseColumn<T1,SO>::ResultType,
                              CompressedVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct AddTrait< CompressedVector<T1,false>, DenseColumn<T2,SO> >
{
   typedef typename AddTrait< CompressedVector<T1,false>,
                              typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< DenseColumn<T1,SO1>, DenseColumn<T2,SO2> >
{
   typedef typename AddTrait< typename DenseColumn<T1,SO1>::ResultType,
                              typename DenseColumn<T2,SO2>::ResultType >::Type  Type;
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
template< typename T1, bool SO, typename T2, size_t N >
struct SubTrait< DenseColumn<T1,SO>, StaticVector<T2,N,false> >
{
   typedef typename SubTrait< typename DenseColumn<T1,SO>::ResultType,
                              StaticVector<T2,N,false> >::Type  Type;
};

template< typename T1, size_t N, typename T2, bool SO >
struct SubTrait< StaticVector<T1,N,false>, DenseColumn<T2,SO> >
{
   typedef typename SubTrait< StaticVector<T1,N,false>,
                              typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< DenseColumn<T1,SO>, DynamicVector<T2,false> >
{
   typedef typename SubTrait< typename DenseColumn<T1,SO>::ResultType,
                              DynamicVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct SubTrait< DynamicVector<T1,false>, DenseColumn<T2,SO> >
{
   typedef typename SubTrait< DynamicVector<T1,false>,
                              typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< DenseColumn<T1,SO>, CompressedVector<T2,false> >
{
   typedef typename SubTrait< typename DenseColumn<T1,SO>::ResultType,
                              CompressedVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct SubTrait< CompressedVector<T1,false>, DenseColumn<T2,SO> >
{
   typedef typename SubTrait< CompressedVector<T1,false>,
                              typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< DenseColumn<T1,SO1>, DenseColumn<T2,SO2> >
{
   typedef typename SubTrait< typename DenseColumn<T1,SO1>::ResultType,
                              typename DenseColumn<T2,SO2>::ResultType >::Type  Type;
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
struct MultTrait< DenseColumn<T1,SO>, T2 >
{
   typedef typename MultTrait< typename DenseColumn<T1,SO>::ResultType, T2 >::Type  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};

template< typename T1, typename T2, bool SO >
struct MultTrait< T1, DenseColumn<T2,SO> >
{
   typedef typename MultTrait< T1, typename DenseColumn<T2,SO>::ResultType >::Type  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T1 );
};

template< typename T1, bool SO, typename T2, size_t N, bool TF >
struct MultTrait< DenseColumn<T1,SO>, StaticVector<T2,N,TF> >
{
   typedef typename MultTrait< typename DenseColumn<T1,SO>::ResultType,
                               StaticVector<T2,N,TF> >::Type  Type;
};

template< typename T1, size_t N, bool TF, typename T2, bool SO >
struct MultTrait< StaticVector<T1,N,TF>, DenseColumn<T2,SO> >
{
   typedef typename MultTrait< StaticVector<T1,N,TF>,
                               typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2, bool TF >
struct MultTrait< DenseColumn<T1,SO>, DynamicVector<T2,TF> >
{
   typedef typename MultTrait< typename DenseColumn<T1,SO>::ResultType,
                               DynamicVector<T2,TF> >::Type  Type;
};

template< typename T1, bool TF, typename T2, bool SO >
struct MultTrait< DynamicVector<T1,TF>, DenseColumn<T2,SO> >
{
   typedef typename MultTrait< DynamicVector<T1,TF>,
                               typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2, bool TF >
struct MultTrait< DenseColumn<T1,SO>, CompressedVector<T2,TF> >
{
   typedef typename MultTrait< typename DenseColumn<T1,SO>::ResultType,
                               CompressedVector<T2,TF> >::Type  Type;
};

template< typename T1, bool TF, typename T2, bool SO >
struct MultTrait< CompressedVector<T1,TF>, DenseColumn<T2,SO> >
{
   typedef typename MultTrait< CompressedVector<T1,TF>,
                               typename DenseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< DenseColumn<T1,SO1>, DenseColumn<T2,SO2> >
{
   typedef typename MultTrait< typename DenseColumn<T1,SO1>::ResultType,
                               typename DenseColumn<T2,SO2>::ResultType >::Type  Type;
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
struct DivTrait< DenseColumn<T1,SO>, T2 >
{
   typedef typename DivTrait< typename DenseColumn<T1,SO>::ResultType, T2 >::Type  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
