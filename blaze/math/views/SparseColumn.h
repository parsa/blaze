//=================================================================================================
/*!
//  \file blaze/math/views/SparseColumn.h
//  \brief Header file for the SparseColumn class template
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

#ifndef _BLAZE_MATH_VIEWS_SPARSECOLUMN_H_
#define _BLAZE_MATH_VIEWS_SPARSECOLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <stdexcept>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/Functions.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
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

//  - MT: specifies the type of the sparse matrix primitive. SparseColumn can be used with any
//        sparse matrix primitive, but does not work with any matrix expression type.
//  - SO: specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the sparse matrix.
//        This template parameter doesn't have to be explicitly defined, but is automatically
//        derived from the first template parameter.
//
// A reference to a sparse column can be conventiently created via the column() function. The row
// can be either used as an alias to grant write access to a specific column of a matrix primitive
// on the left-hand side of an assignment or to grant read-access to a specific column of a matrix
// primitive or expression on the right-hand side of an assignment:

   \code
   blaze::DynamicVector<double,columnVector> x;
   blaze::DynamicMatrix<double,columnMajor> A, B;
   // ... Resizing and initialization

   // Setting the 2nd column of matrix A to x
   column( A, 2UL ) = x;

   // Setting x to the 3rd column of the result of the matrix multiplication
   x = column( A * B, 3UL );
   \endcode

// Inserting/accessing elements in a sparse column can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   using blaze::CompressedMatrix;
   using blaze::SparseColumn;

   typedef CompressedMatrix<double,columnMajor>  MatrixType;
   MatrixType A( 100UL, 10UL );  // Non-initialized 10x100 matrix

   typedef SparseColumn<MatrixType>  ColumnType;
   ColumnType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   // The subscript operator provides access to all possible elements of the sparse column,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse column, the element is inserted into the column.
   col0[42] = 2.0;

   // An alternative for inserting elements into the column is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the column.
   col0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse column is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the column and that the column's capacity is
   // large enough to hold the new element.
   col0.reserve( 10UL );
   col0.append( 51UL, -2.1 );

   // In order to traverse all non-zero elements currently stored in the column, the begin()
   // and end() functions can be used. In the example, all non-zero elements of column are
   // traversed.
   for( ColumnType::Iterator i=a.begin(); i!=a.end(); ++i ) {
      ... = i->value();  // Access to the value of the non-zero element
      ... = i->index();  // Access to the index of the non-zero element
   }
   \endcode

// The following example gives an impression of the use of SparseColumn. All operations (addition,
// subtraction, multiplication, scaling, ...) can be performed on all possible combinations of
// dense and sparse vectors with fitting element types:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::DynamicMatrix;
   using blaze::SparseColumn;

   CompressedVector<double,columnVector> a( 2UL ), b;
   a[1] = 2.0;
   DynamicVector<double,columnVector> c( 2UL, 3UL );

   typedef CompressedMatrix<double,columnMajor>  MatrixType;
   MatrixType A( 2UL, 3UL );  // Non-initialized 2x3 matrix

   typedef SparseColumn<MatrixType>  ColumnType;
   ColumnType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   col0[0] = 0UL;         // Manual initialization of the 0th column of A
   col0[1] = 0UL;
   column( A, 1UL ) = a;  // Dense vector initialization of the 1st column of A
   column( A, 2UL ) = c;  // Sparse vector initialization of the 2nd column of A

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

// It is possible to create a column view on both row-major and column-major matrices. However,
// please note that creating a column view on a matrix stored in row-major fashion can result
// in a considerable performance decrease in comparison to a column view on a column-major matrix
// due to the non-contiguous storage of the non-zero matrix elements. Therefore care has to be
// taken in the choice of the most suitable storage order:

   \code
   // Setup of two row-major matrices
   CompressedMatrix<double,rowMajor> A( 128UL, 128UL );
   CompressedMatrix<double,rowMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th column of the multiplication between A and B ...
   CompressedVector<double,columnVector> x = column( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // A with the 15th column of the row-major matrix B.
   CompressedVector<double,columnVector> x = A * column( B, 15UL );
   \endcode

// Although Blaze performs the resulting matrix/vector multiplication as efficiently as possible
// using a column-major storage order for matrix B would result in a more efficient evaluation.
*/
template< typename MT                                 // Type of the sparse matrix
        , bool SO = IsColumnMajorMatrix<MT>::value >  // Storage order
class SparseColumn : public SparseVector< SparseColumn<MT,SO>, false >
                   , private Expression
{
 private:
   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given sparse matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the sparse column will
       return references and iterators to const. Otherwise \a useConst will be set to 0
       and the sparse column will offer write access to the sparse matrix elements both
       via the subscript operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseColumn<MT,SO>                 This;           //!< Type of this SparseColumn instance.
   typedef typename ColumnTrait<MT>::Type      ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename MT::ElementType            ElementType;    //!< Type of the column elements.
   typedef typename MT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseColumn&                 CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   typedef typename MT::ConstReference  ConstReference;

   //! Reference to a non-constant column value.
   typedef typename SelectType< useConst, ConstReference, typename MT::Reference >::Type  Reference;

   //! Iterator over constant elements.
   typedef typename MT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, typename MT::Iterator >::Type  Iterator;
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
   template< typename VT > inline SparseColumn& operator= ( const DenseVector <VT,false>& rhs );
   template< typename VT > inline SparseColumn& operator= ( const SparseVector<VT,false>& rhs );
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
                              inline ElementType&  insert ( size_t index, const ElementType& value );
                              inline void          erase  ( size_t index );
                              inline Iterator      erase  ( Iterator pos );
                              inline Iterator      find   ( size_t index );
                              inline ConstIterator find   ( size_t index ) const;
                              inline void          reserve( size_t n );
   template< typename Other > inline SparseColumn& scale  ( Other scalar );
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
   MT&          matrix_;  //!< The sparse matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
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
/*!\brief The constructor for SparseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline SparseColumn<MT,SO>::SparseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the column
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
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::Reference SparseColumn<MT,SO>::operator[]( size_t index )
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
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ConstReference SparseColumn<MT,SO>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,col_);
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::Iterator SparseColumn<MT,SO>::begin()
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ConstIterator SparseColumn<MT,SO>::begin() const
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ConstIterator SparseColumn<MT,SO>::cbegin() const
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
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::Iterator SparseColumn<MT,SO>::end()
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ConstIterator SparseColumn<MT,SO>::end() const
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ConstIterator SparseColumn<MT,SO>::cend() const
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
/*!\brief Copy assignment operator for SparseColumn.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::operator=( const SparseColumn& rhs )
{
   using blaze::assign;

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && col_ == rhs.col_ ) )
      return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Column sizes do not match" );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      matrix_.reset  ( col_ );
      matrix_.reserve( col_, tmp.nonZeros() );
      assign( *this, tmp );
   }
   else {
      matrix_.reset  ( col_ );
      matrix_.reserve( col_, rhs.nonZeros() );
      assign( *this, rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different dense vectors.
//
// \param rhs Dense vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::operator=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      matrix_.reset( col_ );
      assign( *this, tmp );
   }
   else {
      matrix_.reset( col_ );
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different sparse vectors.
//
// \param rhs Sparse vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::operator=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &matrix_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      matrix_.reset  ( col_ );
      matrix_.reserve( col_, tmp.nonZeros() );
      assign( *this, tmp );
   }
   else {
      matrix_.reset  ( col_ );
      matrix_.reserve( col_, (~rhs).nonZeros() );
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::operator+=( const Vector<VT,false>& rhs )
{
   using blaze::addAssign;

   if( (~rhs).size() != size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   addAssign( *this, ~rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::operator-=( const Vector<VT,false>& rhs )
{
   using blaze::subAssign;

   if( (~rhs).size() != size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   subAssign( *this, ~rhs );

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::operator*=( const Vector<VT,false>& rhs )
{
   if( (~rhs).size() != size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   typedef typename MultTrait<This,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( MultType );

   const MultType tmp( *this * (~rhs) );
   matrix_.reset( col_ );
   assign( tmp );

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
// This operator can only be used for built-in data types. Additionally, the elements of
// the sparse column must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,SO> >::Type&
   SparseColumn<MT,SO>::operator*=( Other rhs )
{
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
// This operator can only be used for built-in data types. Additionally, the elements of the
// sparse column must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,SO> >::Type&
   SparseColumn<MT,SO>::operator/=( Other rhs )
{
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
        , bool SO >    // Storage order
inline size_t SparseColumn<MT,SO>::size() const
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
        , bool SO >    // Storage order
inline size_t SparseColumn<MT,SO>::capacity() const
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
        , bool SO >    // Storage order
inline size_t SparseColumn<MT,SO>::nonZeros() const
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
        , bool SO >    // Storage order
inline void SparseColumn<MT,SO>::reset()
{
   matrix_.reset( col_ );
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ElementType&
   SparseColumn<MT,SO>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( index, col_, value )->value();
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
        , bool SO >    // Storage order
inline void SparseColumn<MT,SO>::erase( size_t index )
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::Iterator SparseColumn<MT,SO>::erase( Iterator pos )
{
   return matrix_.erase( col_, pos );
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::Iterator SparseColumn<MT,SO>::find( size_t index )
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
        , bool SO >    // Storage order
inline typename SparseColumn<MT,SO>::ConstIterator SparseColumn<MT,SO>::find( size_t index ) const
{
   return matrix_.find( index, col_ );
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
        , bool SO >    // Storage order
void SparseColumn<MT,SO>::reserve( size_t n )
{
   matrix_.reserve( col_, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline SparseColumn<MT,SO>& SparseColumn<MT,SO>::scale( Other scalar )
{
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
template< typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
inline size_t SparseColumn<MT,SO>::extendCapacity() const
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
// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
inline void SparseColumn<MT,SO>::append( size_t index, const ElementType& value, bool check )
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,SO>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool SparseColumn<MT,SO>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,SO>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   size_t nonzeros( 0UL );

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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,SO>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      matrix_.append( element->index(), col_, element->value() );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,SO>::addAssign( const DenseVector<VT,false>& rhs )
{
   typedef typename AddTrait<This,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE       ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( *this + (~rhs) );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,SO>::addAssign( const SparseVector<VT,false>& rhs )
{
   typedef typename AddTrait<This,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE      ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( *this + (~rhs) );
   matrix_.reset  ( col_ );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,SO>::subAssign( const DenseVector<VT,false>& rhs )
{
   typedef typename SubTrait<This,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE       ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( *this - (~rhs) );
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
        , bool SO >      // Storage order
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,SO>::subAssign( const SparseVector<VT,false>& rhs )
{
   typedef typename SubTrait<This,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE      ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( *this - (~rhs) );
   matrix_.reset  ( col_ );
   matrix_.reserve( col_, tmp.nonZeros() );
   assign( tmp );
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SparseColumn for row-major matrices.
// \ingroup views
//
// This specialization of SparseColumn adapts the class template to the requirements of
// row-major matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SparseColumn<MT,false> : public SparseVector< SparseColumn<MT,false>, false >
                             , private Expression
{
 private:
   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given sparse matrix of type
       \a MT is const qualified, \a useConst will be set to 1 and the sparse column will
       return references and iterators to const. Otherwise \a useConst will be set to 0
       and the sparse column will offer write access to the sparse matrix elements both
       via the subscript operator and iterators. */
   enum { useConst = IsConst<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseColumn<MT,false>              This;            //!< Type of this SparseColumn instance.
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

   //**ColumnElement class definition**************************************************************
   /*!\brief Access proxy for a specific element of the sparse column.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnElement
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

      //! Return type for the access to the value of a sparse element.
      typedef typename SelectType< returnConst, ReturnType, ElementType& >::Type  Value;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( matrix )  // The sparse matrix containing the column.
         , row_   ( row    )  // The current row index.
         , column_( column )  // The current column index.
         , pos_   ()          // Iterator to the current sparse element.
      {
         for( ; row_<matrix_.rows(); ++row_ ) {
            pos_ = matrix_.find( row_, column_ );
            if( pos_ != matrix_.end( row_ ) ) break;
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
         : matrix_( matrix )  // The sparse matrix containing the column.
         , row_   ( row    )  // The current row index.
         , column_( column )  // The current column index.
         , pos_   ( pos    )  // Iterator to the current sparse element.
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
         for( ; row_<matrix_.rows(); ++row_ ) {
            pos_ = matrix_.find( row_, column_ );
            if( pos_ != matrix_.end( row_ ) ) break;
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
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const {
         return ( &matrix_ == &rhs.matrix_ ) && ( row_ == rhs.row_ ) && ( column_ == rhs.column_ );
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const {
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
         size_t counter( 0UL );
         for( size_t i=rhs.row_; i<row_; ++i ) {
            if( matrix_.find( i, column_ ) != matrix_.end( i ) )
               ++counter;
         }
         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType&  matrix_;  //!< The sparse matrix containing the column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      /*! \cond BLAZE_INTERNAL */
      template< typename MatrixType2, typename IteratorType2 > friend class ColumnIterator;
      template< typename MT2, bool SO2 > friend class SparseColumn;
      /*! \endcond */
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef ColumnIterator<const MT,typename MT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, ColumnIterator<MT,typename MT::Iterator> >::Type  Iterator;
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
                              inline ElementType&  insert ( size_t index, const ElementType& value );
                              inline void          erase  ( size_t index );
                              inline Iterator      erase  ( Iterator pos );
                              inline Iterator      find   ( size_t index );
                              inline ConstIterator find   ( size_t index ) const;
                              inline void          reserve( size_t n );
   template< typename Other > inline SparseColumn& scale  ( Other scalar );
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
   MT&          matrix_;  //!< The sparse matrix containing the column.
   const size_t col_;     //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
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
/*!\brief The constructor for SparseColumn.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseColumn<MT,false>::SparseColumn( MT& matrix, size_t index )
   : matrix_( matrix )  // The sparse matrix containing the column
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
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false>::Reference SparseColumn<MT,false>::operator[]( size_t index )
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
template< typename MT >  // Type of the sparse matrix
inline typename SparseColumn<MT,false>::ConstReference SparseColumn<MT,false>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,col_);
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
inline typename SparseColumn<MT,false>::Iterator SparseColumn<MT,false>::begin()
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
inline typename SparseColumn<MT,false>::ConstIterator SparseColumn<MT,false>::begin() const
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
inline typename SparseColumn<MT,false>::ConstIterator SparseColumn<MT,false>::cbegin() const
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
inline typename SparseColumn<MT,false>::Iterator SparseColumn<MT,false>::end()
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
inline typename SparseColumn<MT,false>::ConstIterator SparseColumn<MT,false>::end() const
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
inline typename SparseColumn<MT,false>::ConstIterator SparseColumn<MT,false>::cend() const
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
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
inline SparseColumn<MT,false>& SparseColumn<MT,false>::operator=( const SparseColumn& rhs )
{
   using blaze::assign;

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && col_ == rhs.col_ ) )
      return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Column sizes do not match" );

   if( rhs.canAlias( &matrix_ ) ) {
      const ResultType tmp( rhs );
      assign( *this, tmp );
   }
   else {
      assign( *this, rhs );
   }

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
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false>& SparseColumn<MT,false>::operator=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   const typename VT::CompositeType tmp( ~rhs );
   assign( *this, tmp );

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false>& SparseColumn<MT,false>::operator+=( const Vector<VT,false>& rhs )
{
   using blaze::addAssign;

   if( (~rhs).size() != size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   addAssign( *this, ~rhs );

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false>& SparseColumn<MT,false>::operator-=( const Vector<VT,false>& rhs )
{
   using blaze::subAssign;

   if( (~rhs).size() != size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   subAssign( *this, ~rhs );

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
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side vector
inline SparseColumn<MT,false>& SparseColumn<MT,false>::operator*=( const Vector<VT,false>& rhs )
{
   if( (~rhs).size() != size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   typedef typename MultTrait<This,typename VT::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( MultType );

   const MultType tmp( *this * (~rhs) );
   assign( tmp );

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
// This operator can only be used for built-in data types. Additionally, the elements of
// the sparse column must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,false> >::Type&
   SparseColumn<MT,false>::operator*=( Other rhs )
{
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
// This operator can only be used for built-in data types. Additionally, the elements of the
// sparse column must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseColumn<MT,false> >::Type&
   SparseColumn<MT,false>::operator/=( Other rhs )
{
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
inline size_t SparseColumn<MT,false>::size() const
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
inline size_t SparseColumn<MT,false>::capacity() const
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
inline size_t SparseColumn<MT,false>::nonZeros() const
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
inline void SparseColumn<MT,false>::reset()
{
   for( size_t i=0UL; i<size(); ++i ) {
      matrix_.erase( i, col_ );
   }
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
inline typename SparseColumn<MT,false>::ElementType&
   SparseColumn<MT,false>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( index, col_, value )->value();
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
inline void SparseColumn<MT,false>::erase( size_t index )
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
inline typename SparseColumn<MT,false>::Iterator SparseColumn<MT,false>::erase( Iterator pos )
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
inline typename SparseColumn<MT,false>::Iterator SparseColumn<MT,false>::find( size_t index )
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
inline typename SparseColumn<MT,false>::ConstIterator SparseColumn<MT,false>::find( size_t index ) const
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
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT >  // Type of the sparse matrix
void SparseColumn<MT,false>::reserve( size_t n )
{
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
*/
template< typename MT >     // Type of the sparse matrix
template< typename Other >  // Data type of the scalar value
inline SparseColumn<MT,false>& SparseColumn<MT,false>::scale( Other scalar )
{
   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
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
// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT >  // Type of the sparse matrix
inline void SparseColumn<MT,false>::append( size_t index, const ElementType& value, bool /*check*/ )
{
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
inline bool SparseColumn<MT,false>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &matrix_ ) == static_cast<const void*>( alias );
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
inline bool SparseColumn<MT,false>::isAliased( const Other* alias ) const
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
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side dense vector
inline void SparseColumn<MT,false>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      const typename MT::Iterator pos = matrix_.insert( i, col_, (~rhs)[i] );
      if( isDefault( pos->value() ) ) {
         matrix_.erase( i, pos );
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
template< typename MT >  // Type of the sparse matrix
template< typename VT >  // Type of the right-hand side sparse vector
inline void SparseColumn<MT,false>::assign( const SparseVector<VT,false>& rhs )
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
inline void SparseColumn<MT,false>::addAssign( const Vector<VT,false>& rhs )
{
   typedef typename AddTrait<This,typename VT::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( *this + (~rhs) );
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
inline void SparseColumn<MT,false>::subAssign( const Vector<VT,false>& rhs )
{
   typedef typename SubTrait<This,typename VT::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION     ( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( *this - (~rhs) );
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
template< typename MT, bool SO >
inline void reset( SparseColumn<MT,SO>& column );

template< typename MT, bool SO >
inline void clear( SparseColumn<MT,SO>& column );

template< typename MT, bool SO >
inline bool isnan( const SparseColumn<MT,SO>& column );

template< typename MT, bool SO >
inline bool isDefault( const SparseColumn<MT,SO>& column );
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
        , bool SO >    // Storage order
inline void reset( SparseColumn<MT,SO>& column )
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
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline void clear( SparseColumn<MT,SO>& column )
{
   column.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse column for not-a-number elements.
// \ingroup sparse_column
//
// \param column The sparse column to be checked for not-a-number elements.
// \return \a true if at least one element of the column is not-a-number, \a false otherwise.
//
// This function checks the sparse column for not-a-number (NaN) elements. If at least one element
// of the column is not-a-number, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isnan( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline bool isnan( const SparseColumn<MT,SO>& column )
{
   typedef typename SparseColumn<MT,SO>::ConstIterator  ConstIterator;

   const ConstIterator end( column.end() );
   for( ConstIterator element=column.begin(); element!=end; ++element ) {
      if( isnan( element->value() ) ) return true;
   }
   return false;
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
        , bool SO >    // Storage order
inline bool isDefault( const SparseColumn<MT,SO>& column )
{
   typedef typename SparseColumn<MT,SO>::ConstIterator  ConstIterator;

   const ConstIterator end( column.end() );
   for( ConstIterator element=column.begin(); element!=end; ++element )
      if( !isDefault( element->value() ) ) return false;
   return true;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given sparse matrix.
// \ingroup views
//
// \param sm The sparse matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given sparse matrix.

   \code
   typedef blaze::CompressedMatrix<double,columnMajor>  Matrix;

   Matrix A;
   // ... Resizing and initialization
   blaze::SparseColumn<Matrix> = column( A, 3UL );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline SparseColumn<MT> column( SparseMatrix<MT,SO>& sm, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return SparseColumn<MT>( ~sm, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant sparse matrix.
// \ingroup views
//
// \param sm The constant sparse matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given sparse matrix.

   \code
   typedef blaze::CompressedMatrix<double,columnMajor>  Matrix;

   const Matrix A;
   // ... Resizing and initialization
   blaze::SparseColumn<const Matrix> = column( A, 3UL );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline SparseColumn<const MT> column( const SparseMatrix<MT,SO>& sm, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return SparseColumn<const MT>( ~sm, index );
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
struct AddTrait< SparseColumn<T1,SO>, StaticVector<T2,N,false> >
{
   typedef typename AddTrait< typename SparseColumn<T1,SO>::ResultType,
                              StaticVector<T2,N,false> >::Type  Type;
};

template< typename T1, size_t N, typename T2, bool SO >
struct AddTrait< StaticVector<T1,N,false>, SparseColumn<T2,SO> >
{
   typedef typename AddTrait< StaticVector<T1,N,false>,
                              typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< SparseColumn<T1,SO>, DynamicVector<T2,false> >
{
   typedef typename AddTrait< typename SparseColumn<T1,SO>::ResultType,
                              DynamicVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct AddTrait< DynamicVector<T1,false>, SparseColumn<T2,SO> >
{
   typedef typename AddTrait< DynamicVector<T1,false>,
                              typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< SparseColumn<T1,SO1>, DenseColumn<T2,SO2> >
{
   typedef typename AddTrait< typename SparseColumn<T1,SO1>::ResultType,
                              typename DenseColumn <T2,SO2>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< DenseColumn<T1,SO1>, SparseColumn<T2,SO2> >
{
   typedef typename AddTrait< typename DenseColumn <T1,SO1>::ResultType,
                              typename SparseColumn<T2,SO2>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< SparseColumn<T1,SO>, CompressedVector<T2,false> >
{
   typedef typename AddTrait< typename SparseColumn<T1,SO>::ResultType,
                              CompressedVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct AddTrait< CompressedVector<T1,false>, SparseColumn<T2,SO> >
{
   typedef typename AddTrait< CompressedVector<T1,false>,
                              typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< SparseColumn<T1,SO1>, SparseColumn<T2,SO2> >
{
   typedef typename AddTrait< typename SparseColumn<T1,SO1>::ResultType,
                              typename SparseColumn<T2,SO2>::ResultType >::Type  Type;
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
struct SubTrait< SparseColumn<T1,SO>, StaticVector<T2,N,false> >
{
   typedef typename SubTrait< typename SparseColumn<T1,SO>::ResultType,
                              StaticVector<T2,N,false> >::Type  Type;
};

template< typename T1, size_t N, typename T2, bool SO >
struct SubTrait< StaticVector<T1,N,false>, SparseColumn<T2,SO> >
{
   typedef typename SubTrait< StaticVector<T1,N,false>,
                              typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< SparseColumn<T1,SO>, DynamicVector<T2,false> >
{
   typedef typename SubTrait< typename SparseColumn<T1,SO>::ResultType,
                              DynamicVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct SubTrait< DynamicVector<T1,false>, SparseColumn<T2,SO> >
{
   typedef typename SubTrait< DynamicVector<T1,false>,
                              typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< SparseColumn<T1,SO1>, DenseColumn<T2,SO2> >
{
   typedef typename SubTrait< typename SparseColumn<T1,SO1>::ResultType,
                              typename DenseColumn <T2,SO2>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< DenseColumn<T1,SO1>, SparseColumn<T2,SO2> >
{
   typedef typename SubTrait< typename DenseColumn <T1,SO1>::ResultType,
                              typename SparseColumn<T2,SO2>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< SparseColumn<T1,SO>, CompressedVector<T2,false> >
{
   typedef typename SubTrait< typename SparseColumn<T1,SO>::ResultType,
                              CompressedVector<T2,false> >::Type  Type;
};

template< typename T1, typename T2, bool SO >
struct SubTrait< CompressedVector<T1,false>, SparseColumn<T2,SO> >
{
   typedef typename SubTrait< CompressedVector<T1,false>,
                              typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< SparseColumn<T1,SO1>, SparseColumn<T2,SO2> >
{
   typedef typename SubTrait< typename SparseColumn<T1,SO1>::ResultType,
                              typename SparseColumn<T2,SO2>::ResultType >::Type  Type;
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
struct MultTrait< SparseColumn<T1,SO>, T2 >
{
   typedef typename MultTrait< typename SparseColumn<T1,SO>::ResultType, T2 >::Type  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};

template< typename T1, typename T2, bool SO >
struct MultTrait< T1, SparseColumn<T2,SO> >
{
   typedef typename MultTrait< T1, typename SparseColumn<T2,SO>::ResultType >::Type  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T1 );
};

template< typename T1, bool SO, typename T2, size_t N, bool TF >
struct MultTrait< SparseColumn<T1,SO>, StaticVector<T2,N,TF> >
{
   typedef typename MultTrait< typename SparseColumn<T1,SO>::ResultType,
                               StaticVector<T2,N,TF> >::Type  Type;
};

template< typename T1, size_t N, bool TF, typename T2, bool SO >
struct MultTrait< StaticVector<T1,N,TF>, SparseColumn<T2,SO> >
{
   typedef typename MultTrait< StaticVector<T1,N,TF>,
                               typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2, bool TF >
struct MultTrait< SparseColumn<T1,SO>, DynamicVector<T2,TF> >
{
   typedef typename MultTrait< typename SparseColumn<T1,SO>::ResultType,
                               DynamicVector<T2,TF> >::Type  Type;
};

template< typename T1, bool TF, typename T2, bool SO >
struct MultTrait< DynamicVector<T1,TF>, SparseColumn<T2,SO> >
{
   typedef typename MultTrait< DynamicVector<T1,TF>,
                               typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< SparseColumn<T1,SO1>, DenseColumn<T2,SO2> >
{
   typedef typename MultTrait< typename SparseColumn<T1,SO1>::ResultType,
                               typename DenseColumn <T2,SO2>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< DenseColumn<T1,SO1>, SparseColumn<T2,SO2> >
{
   typedef typename MultTrait< typename DenseColumn <T1,SO1>::ResultType,
                               typename SparseColumn<T2,SO2>::ResultType >::Type  Type;
};

template< typename T1, bool SO, typename T2, bool TF >
struct MultTrait< SparseColumn<T1,SO>, CompressedVector<T2,TF> >
{
   typedef typename MultTrait< typename SparseColumn<T1,SO>::ResultType,
                               CompressedVector<T2,TF> >::Type  Type;
};

template< typename T1, bool TF, typename T2, bool SO >
struct MultTrait< CompressedVector<T1,TF>, SparseColumn<T2,SO> >
{
   typedef typename MultTrait< CompressedVector<T1,TF>,
                               typename SparseColumn<T2,SO>::ResultType >::Type  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< SparseColumn<T1,SO1>, SparseColumn<T2,SO2> >
{
   typedef typename MultTrait< typename SparseColumn<T1,SO1>::ResultType,
                               typename SparseColumn<T2,SO2>::ResultType >::Type  Type;
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
struct DivTrait< SparseColumn<T1,SO>, T2 >
{
   typedef typename DivTrait< typename SparseColumn<T1,SO>::ResultType, T2 >::Type  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
