//=================================================================================================
/*!
//  \file blaze/math/views/DenseSubvector.h
//  \brief Header file for the DenseSubvector class template
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

#ifndef _BLAZE_MATH_VIEWS_DENSESUBVECTOR_H_
#define _BLAZE_MATH_VIEWS_DENSESUBVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/CrossExpr.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/Forward.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsCrossExpr.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/math/typetraits/IsMatVecMultExpr.h>
#include <blaze/math/typetraits/IsTVecMatMultExpr.h>
#include <blaze/math/typetraits/IsVecAbsExpr.h>
#include <blaze/math/typetraits/IsVecEvalExpr.h>
#include <blaze/math/typetraits/IsVecScalarDivExpr.h>
#include <blaze/math/typetraits/IsVecScalarMultExpr.h>
#include <blaze/math/typetraits/IsVecTransExpr.h>
#include <blaze/math/typetraits/IsVecVecAddExpr.h>
#include <blaze/math/typetraits/IsVecVecMultExpr.h>
#include <blaze/math/typetraits/IsVecVecSubExpr.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/CacheSize.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_subvector Dense Subvector
// \ingroup views
*/
/*!\brief View on a specific subvector of a dense vector.
// \ingroup dense_subvector
//
// The DenseSubvector template represents a view on a specific subvector of a dense vector
// primitive. The type of the dense vector is specified via the first template parameter:

   \code
   template< typename VT, bool TF >
   class DenseSubvector;
   \endcode

//  - VT: specifies the type of the dense vector primitive. DenseSubvector can be used with every
//        dense vector primitive or view, but does not work with any vector expression type.
//  - TF: specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//        vector (\a blaze::columnVector). This template parameter doesn't have to be explicitly
//        defined, but is automatically derived from the first template parameter.
//
//
// \n \section dense_subvector_setup Setup of Dense Subvectors
//
// A view on a dense subvector can be created very conveniently via the \c subvector() function.
// This view can be treated as any other dense vector, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The view can also be used on both
// sides of an assignment: The subvector can either be used as an alias to grant write access to
// a specific subvector of a dense vector primitive on the left-hand side of an assignment or
// to grant read-access to a specific subvector of a dense vector primitive or expression on
// the right-hand side of an assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>      DenseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   DenseMatrixType  A;
   // ... Resizing and initialization

   // Setting the first ten elements of x to the 2nd row of matrix A
   blaze::DenseSubvector<DenseVectorType> sv = subvector( x, 0UL, 10UL );
   sv = row( A, 2UL );

   // Setting the second ten elements of x to y
   subvector( x, 10UL, 10UL ) = y;

   // Setting the 3rd row of A to a subvector of x
   row( A, 3UL ) = subvector( x, 3UL, 10UL );

   // Setting x to a subvector of the result of the addition between y and the 1st row of A
   x = subvector( y + row( A, 1UL ), 2UL, 5UL )
   \endcode

// \n \section dense_subvector_element_access Element access
//
// A dense subvector can be used like any other dense vector. For instance, the elements of the
// dense subvector can be directly accessed with the subscript operator.

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>  VectorType;
   VectorType v;
   // ... Resizing and initialization

   // Creating an 8-dimensional subvector, starting from index 4
   blaze::DenseSubvector<VectorType> sv = subvector( v, 4UL, 8UL );

   // Setting the 1st element of the subvector, which corresponds to
   // the element at index 5 in vector v
   sv[1] = 2.0;
   \endcode

// The numbering of the subvector elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the specified size of the subvector. Alternatively, the elements of a subvector can
// be traversed via iterators. Just as with vectors, in case of non-const subvectors, \c begin()
// and \c end() return an Iterator, which allows a manipulation of the non-zero values, in case
// of constant subvectors a ConstIterator is returned:

   \code
   typedef blaze::DynamicVector<int,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>           SubvectorType;

   VectorType v( 256UL );
   // ... Resizing and initialization

   // Creating a reference to a specific subvector of vector v
   SubvectorType sv = subvector( v, 16UL, 64UL );

   for( SubvectorType::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
      *it = ...;  // OK: Write access to the dense subvector value.
      ... = *it;  // OK: Read access to the dense subvector value.
   }

   for( SubvectorType::ConstIterator it=sv.begin(); it!=sv.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense subvector value.
   }
   \endcode

// \n \section dense_subvector_common_operations Common Operations
//
// The current number of subvector elements can be obtained via the \c size() function, the
// current capacity via the \c capacity() function, and the number of non-zero elements via
// the \c nonZeros() function. However, since subvector are views on a specific subvector of
// a vector, several operations are not possible on views, such as resizing and swapping:

   \code
   typedef blaze::DynamicVector<int,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>           SubvectorType;

   VectorType v( 42UL );
   // ... Resizing and initialization

   // Creating a view on the range [5..15] of vector v
   SubvectorType sv = subvector( v, 5UL, 10UL );

   sv.size();          // Returns the number of elements in the subvector
   sv.capacity();      // Returns the capacity of the subvector
   sv.nonZeros();      // Returns the number of non-zero elements contained in the subvector

   sv.resize( 84UL );  // Compilation error: Cannot resize a subvector of a vector

   SubvectorType sv2 = subvector( v, 15UL, 10UL );
   swap( sv, sv2 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section sparse_subvector_arithmetic_operations Arithmetic Operations
//
// The following example gives an impression of the use of DenseSubvector within arithmetic
// operations. All operations (addition, subtraction, multiplication, scaling, ...) can be
// performed on all possible combinations of dense and sparse vectors with fitting element
// types:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   DenseVectorType d1, d2, d3;
   SparseVectorType s1, s2;

   // ... Resizing and initialization

   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrixType;
   DenseMatrixType A;

   typedef blaze::DenseSubvector<DenseVectorType>  SubvectorType;
   SubvectorType sv( subvector( d1, 0UL, 10UL ) );  // View on the range [0..9] of vector d1

   sv = d2;                           // Dense vector initialization of the range [0..9]
   subvector( d1, 10UL, 10UL ) = s1;  // Sparse vector initialization of the range [10..19]

   d3 = sv + d2;                           // Dense vector/dense vector addition
   s2 = s1 + subvector( d1, 10UL, 10UL );  // Sparse vector/dense vector addition
   d2 = sv * subvector( d1, 20UL, 10UL );  // Component-wise vector multiplication

   subvector( d1, 3UL, 4UL ) *= 2.0;      // In-place scaling of the range [3..6]
   d2 = subvector( d1, 7UL, 3UL ) * 2.0;  // Scaling of the range [7..9]
   d2 = 2.0 * subvector( d1, 7UL, 3UL );  // Scaling of the range [7..9]

   subvector( d1, 0UL , 10UL ) += d2;  // Addition assignment
   subvector( d1, 10UL, 10UL ) -= s2;  // Subtraction assignment
   subvector( d1, 20UL, 10UL ) *= sv;  // Multiplication assignment

   double scalar = subvector( d1, 5UL, 10UL ) * trans( s1 );  // Scalar/dot/inner product between two vectors

   A = trans( s1 ) * subvector( d1, 4UL, 16UL );  // Outer product between two vectors
   \endcode

// \n \section dense_subvector_on_dense_subvector Subvectors on Subvectors
//
// It is also possible to create a subvector view on another subvector. In this context it is
// important to remember that the type returned by the \c subvector() function is the same type
// as the type of the given subvector, since the view on a subvector is just another view on the
// underlying dense vector:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>              SubvectorType;

   VectorType d1;

   // ... Resizing and initialization

   // Creating a subvector view on the dense vector d1
   SubvectorType sv1 = subvector( d1, 5UL, 10UL );

   // Creating a subvector view on the dense subvector sv1
   SubvectorType sv2 = subvector( sv1, 1UL, 5UL );
   \endcode
*/
template< typename VT                               // Type of the dense vector
        , bool TF = IsTransposeVector<VT>::value >  // Transpose flag
class DenseSubvector : public DenseVector< DenseSubvector<VT,TF>, TF >
                     , private Expression
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, VT, VT& >::Type  Operand;

   //! Intrinsic trait for the vector element type.
   typedef IntrinsicTrait<typename VT::ElementType>  IT;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the non-const reference and iterator types.
   /*! The \a useConst compile time constant expression represents a compilation switch for
       the non-const reference and iterator types. In case the given dense vector of type
       \a VT is const qualified, \a useConst will be set to 1 and the subvector will return
       references and iterators to const. Otherwise \a useConst will be set to 0 and the
       subvector will offer write access to the dense vector elements both via the subscript
       operator and iterators. */
   enum { useConst = IsConst<VT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<VT,TF>               This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<VT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename VT::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the subvector elements.
   typedef typename VT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubvector&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   typedef typename VT::ConstReference  ConstReference;

   //! Reference to a non-constant subvector value.
   typedef typename SelectType< useConst, ConstReference, typename VT::Reference >::Type  Reference;

   //! Pointer to a constant subvector value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a constant subvector value.
   typedef typename SelectType< useConst, ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename VT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename SelectType< useConst, ConstIterator, typename VT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT::vectorizable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubvector( VT& vector, size_t index, size_t n );
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
                            inline DenseSubvector& operator= ( const ElementType& rhs );
                            inline DenseSubvector& operator= ( const DenseSubvector& rhs );
   template< typename VT2 > inline DenseSubvector& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator*=( const Vector<VT2,TF>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubvector >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DenseSubvector >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          size() const;
                              inline size_t          capacity() const;
                              inline size_t          nonZeros() const;
                              inline void            reset();
   template< typename Other > inline DenseSubvector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedAssign {
      enum { value = vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedAddAssign {
      enum { value = vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedSubAssign {
      enum { value = vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedMultAssign {
      enum { value = vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
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

   inline IntrinsicType load  ( size_t index ) const;
   inline IntrinsicType loadu ( size_t index ) const;
   inline void          store ( size_t index, const IntrinsicType& value );
   inline void          storeu( size_t index, const IntrinsicType& value );
   inline void          stream( size_t index, const IntrinsicType& value );

   template< typename VT2 >
   inline typename DisableIf< VectorizedAssign<VT2> >::Type
      assign( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 >
   inline typename EnableIf< VectorizedAssign<VT2> >::Type
      assign( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 > inline void assign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline typename DisableIf< VectorizedAddAssign<VT2> >::Type
      addAssign( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 >
   inline typename EnableIf< VectorizedAddAssign<VT2> >::Type
      addAssign ( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 > inline void addAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline typename DisableIf< VectorizedSubAssign<VT2> >::Type
      subAssign ( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 >
   inline typename EnableIf< VectorizedSubAssign<VT2> >::Type
      subAssign( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 > inline void subAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline typename DisableIf< VectorizedMultAssign<VT2> >::Type
      multAssign( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 >
   inline typename EnableIf< VectorizedMultAssign<VT2> >::Type
      multAssign( const DenseVector <VT2,TF>& rhs );

   template< typename VT2 > inline void multAssign( const SparseVector<VT2,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      vector_;   //!< The dense vector containing the subvector.
   const size_t offset_;   //!< The offset of the subvector within the dense vector.
   const size_t size_;     //!< The size of the subvector.
   const size_t rest_;     //!< The number of remaining elements in an unaligned intrinsic operation.
   const size_t final_;    //!< The final index for unaligned intrinsic operations.
                           /*!< In case the subvector is not fully aligned and the subvector is
                                involved in a vectorized operation, the final index indicates at
                                which index a special treatment for the remaining elements is
                                required. */
   const bool   aligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the subvector is fully aligned.
                                In case the subvector is fully aligned, no special handling has to
                                be used for the last elements of the subvector in a vectorized
                                operation. In order to be aligned, the following conditions must
                                hold for the subvector:
                                 - The first element of the subvector must be aligned
                                 - The subvector must be at the end of the given vector or
                                 - The size of the subvector must be a multiple of the number of
                                   values per intrinsic element. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename VT2, bool TF2 >
   friend DenseSubvector<VT2,TF2> subvector( const DenseSubvector<VT2,TF2>& dv, size_t index, size_t size );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
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
/*!\brief The constructor for DenseSubvector.
//
// \param vector The dense vector containing the subvector.
// \param index The first index of the subvector in the given vector.
// \param n The size of the subvector.
// \exception std::invalid_argument Invalid subvector specification.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,TF>::DenseSubvector( VT& vector, size_t index, size_t n )
   : vector_ ( vector       )  // The vector containing the subvector
   , offset_ ( index        )  // The offset of the subvector within the sparse vector
   , size_   ( n            )  // The size of the subvector
   , rest_   ( n % IT::size )  // The number of remaining elements in an unaligned intrinsic operation
   , final_  ( n - rest_    )  // The final index for unaligned intrinsic operations
   , aligned_( ( index % IT::size == 0UL ) &&
               ( index + n == vector.size() || n % IT::size == 0UL ) )
{
   if( index + n > vector.size() )
      throw std::invalid_argument( "Invalid subvector specification" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::Reference
   DenseSubvector<VT,TF>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return vector_[offset_+index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::ConstReference
   DenseSubvector<VT,TF>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return const_cast<const VT&>( vector_ )[offset_+index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::Pointer DenseSubvector<VT,TF>::data()
{
   return vector_.data() + offset_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::ConstPointer DenseSubvector<VT,TF>::data() const
{
   return vector_.data() + offset_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::Iterator DenseSubvector<VT,TF>::begin()
{
   return vector_.begin() + offset_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::ConstIterator DenseSubvector<VT,TF>::begin() const
{
   return vector_.begin() + offset_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::ConstIterator DenseSubvector<VT,TF>::cbegin() const
{
   return vector_.begin() + offset_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::Iterator DenseSubvector<VT,TF>::end()
{
   return vector_.begin() + offset_ + size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::ConstIterator DenseSubvector<VT,TF>::end() const
{
   return vector_.begin() + offset_ + size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::ConstIterator DenseSubvector<VT,TF>::cend() const
{
   return vector_.begin() + offset_ + size_;
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Homogenous assignment to all subvector elements.
//
// \param rhs Scalar value to be assigned to all subvector elements.
// \return Reference to the assigned subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::operator=( const ElementType& rhs )
{
   const size_t iend( offset_ + size_ );

   for( size_t i=offset_; i<iend; ++i )
      vector_[i] = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DenseSubvector.
//
// \param rhs Dense subvector to be copied.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Subvector sizes do not match.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::operator=( const DenseSubvector& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( &rhs == this || ( &vector_ == &rhs.vector_ && offset_ == rhs.offset_ ) )
      return *this;

   if( size() != rhs.size() )
      throw std::invalid_argument( "Subvector sizes do not match" );

   if( rhs.canAlias( &vector_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      assign( *this, tmp );
   }
   else {
      assign( *this, rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::operator=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &vector_ ) ) {
      const typename VT::ResultType tmp( ~rhs );
      assign( *this, tmp );
   }
   else {
      if( IsSparseVector<VT2>::value )
         reset();
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::operator+=( const Vector<VT2,TF>& rhs )
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &vector_ ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::operator-=( const Vector<VT2,TF>& rhs )
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &vector_ ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::operator*=( const Vector<VT2,TF>& rhs )
{
   using blaze::multAssign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT::ResultType );

   if( size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~rhs).canAlias( &vector_ ) || IsSparseVector<VT2>::value ) {
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
/*!\brief Multiplication assignment operator for the multiplication between a subvector and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the assigned subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubvector<VT,TF> >::Type&
   DenseSubvector<VT,TF>::operator*=( Other rhs )
{
   return operator=( (*this) * rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a subvector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the assigned subvector.
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubvector<VT,TF> >::Type&
   DenseSubvector<VT,TF>::operator/=( Other rhs )
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
/*!\brief Returns the current size/dimension of the dense subvector.
//
// \return The size of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,TF>::size() const
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the dense subvector.
//
// \return The capacity of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,TF>::capacity() const
{
   return vector_.capacity() - offset_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the subvector.
//
// \return The number of non-zero elements in the subvector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,TF>::nonZeros() const
{
   size_t nonzeros( 0 );

   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i ) {
      if( !isDefault( vector_[i] ) )
         ++nonzeros;
   }

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,TF>::reset()
{
   using blaze::reset;

   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i )
      reset( vector_[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the dense subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the scalar value
inline DenseSubvector<VT,TF>& DenseSubvector<VT,TF>::scale( const Other& scalar )
{
   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i )
      vector_[i] *= scalar;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the dense subvector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubvector<VT,TF>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( &vector_ ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense subvector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubvector<VT,TF>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( &vector_ ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of an intrinsic element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::IntrinsicType
   DenseSubvector<VT,TF>::load( size_t index ) const
{
   return loadu( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of an intrinsic element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense
// subvector. The index must be smaller than the number of subvector elements and it must be
// a multiple of the number of values inside the intrinsic element. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,TF>::IntrinsicType
   DenseSubvector<VT,TF>::loadu( size_t index ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()         , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL, "Invalid subvector access index" );

   if( aligned_ || index != final_ ) {
      return vector_.loadu( offset_+index );
   }
   else {
      IntrinsicType value;
      for( size_t j=0UL; j<rest_; ++j )
         value[j] = vector_[offset_+index+j];
      return value;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of an intrinsic element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store a specific intrinsic element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,TF>::store( size_t index, const IntrinsicType& value )
{
   storeu( index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of an intrinsic element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store a specific intrinsic element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,TF>::storeu( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()         , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL, "Invalid subvector access index" );

   if( aligned_ || index != final_ ) {
      vector_.storeu( offset_+index, value );
   }
   else {
      for( size_t j=0UL; j<rest_; ++j )
         vector_[offset_+index+j] = value[j];
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of an intrinsic element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific intrinsic element of the
// dense subvector. The index must be smaller than the number of subvector elements and it
// must be a multiple of the number of values inside the intrinsic element. This function
// must \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result in
// erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,TF>::stream( size_t index, const IntrinsicType& value )
{
   storeu( index, value );
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedAssign<VT2> >::Type
   DenseSubvector<VT,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      vector_[i+offset_    ] = (~rhs)[i    ];
      vector_[i+offset_+1UL] = (~rhs)[i+1UL];
   }
   if( iend < size() ) {
      vector_[iend+offset_] = (~rhs)[iend];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedAssign<VT2> >::Type
   DenseSubvector<VT,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   if( aligned_ && ( size_ > ( cacheSize/( sizeof(ElementType) * 3UL ) ) ) && !(~rhs).isAliased( &vector_ ) )
   {
      for( size_t i=0UL; i<size(); i+=IT::size ) {
         vector_.stream( offset_+i, (~rhs).load(i) );
      }
   }
   else
   {
      const size_t iend( size_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

      for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
         vector_.storeu( offset_+i             , (~rhs).load(i             ) );
         vector_.storeu( offset_+i+IT::size    , (~rhs).load(i+IT::size    ) );
         vector_.storeu( offset_+i+IT::size*2UL, (~rhs).load(i+IT::size*2UL) );
         vector_.storeu( offset_+i+IT::size*3UL, (~rhs).load(i+IT::size*3UL) );
      }
      for( size_t i=iend; i<size_; i+=IT::size ) {
         storeu( i, (~rhs).load(i) );
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,TF>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[element->index()+offset_] = element->value();
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedAddAssign<VT2> >::Type
   DenseSubvector<VT,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      vector_[i+offset_    ] += (~rhs)[i    ];
      vector_[i+offset_+1UL] += (~rhs)[i+1UL];
   }
   if( iend < size() ) {
      vector_[iend+offset_] += (~rhs)[iend];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedAddAssign<VT2> >::Type
   DenseSubvector<VT,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t iend( size_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

   for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
      vector_.storeu( offset_+i             , load(i             ) + (~rhs).load(i             ) );
      vector_.storeu( offset_+i+IT::size    , load(i+IT::size    ) + (~rhs).load(i+IT::size    ) );
      vector_.storeu( offset_+i+IT::size*2UL, load(i+IT::size*2UL) + (~rhs).load(i+IT::size*2UL) );
      vector_.storeu( offset_+i+IT::size*3UL, load(i+IT::size*3UL) + (~rhs).load(i+IT::size*3UL) );
   }
   for( size_t i=iend; i<size_; i+=IT::size ) {
      storeu( i, load(i) + (~rhs).load(i) );
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,TF>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[element->index()+offset_] += element->value();
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedSubAssign<VT2> >::Type
   DenseSubvector<VT,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      vector_[i+offset_    ] -= (~rhs)[i    ];
      vector_[i+offset_+1UL] -= (~rhs)[i+1UL];
   }
   if( iend < size() ) {
      vector_[iend+offset_] -= (~rhs)[iend];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedSubAssign<VT2> >::Type
   DenseSubvector<VT,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t iend( size_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

   for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
      vector_.storeu( offset_+i             , load(i             ) - (~rhs).load(i             ) );
      vector_.storeu( offset_+i+IT::size    , load(i+IT::size    ) - (~rhs).load(i+IT::size    ) );
      vector_.storeu( offset_+i+IT::size*2UL, load(i+IT::size*2UL) - (~rhs).load(i+IT::size*2UL) );
      vector_.storeu( offset_+i+IT::size*3UL, load(i+IT::size*3UL) - (~rhs).load(i+IT::size*3UL) );
   }
   for( size_t i=iend; i<size_; i+=IT::size ) {
      storeu( i, load(i) - (~rhs).load(i) );
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,TF>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[element->index()+offset_] -= element->value();
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedMultAssign<VT2> >::Type
   DenseSubvector<VT,TF>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t iend( size() & size_t(-2) );
   for( size_t i=0UL; i<iend; i+=2UL ) {
      vector_[i+offset_    ] *= (~rhs)[i    ];
      vector_[i+offset_+1UL] *= (~rhs)[i+1UL];
   }
   if( iend < size() ) {
      vector_[iend+offset_] *= (~rhs)[iend];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,TF>::BLAZE_TEMPLATE VectorizedMultAssign<VT2> >::Type
   DenseSubvector<VT,TF>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   const size_t iend( size_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

   for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
      vector_.storeu( offset_+i             , load(i             ) * (~rhs).load(i             ) );
      vector_.storeu( offset_+i+IT::size    , load(i+IT::size    ) * (~rhs).load(i+IT::size    ) );
      vector_.storeu( offset_+i+IT::size*2UL, load(i+IT::size*2UL) * (~rhs).load(i+IT::size*2UL) );
      vector_.storeu( offset_+i+IT::size*3UL, load(i+IT::size*3UL) * (~rhs).load(i+IT::size*3UL) );
   }
   for( size_t i=iend; i<size_; i+=IT::size ) {
      storeu( i, load(i) * (~rhs).load(i) );
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,TF>::multAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( *this );

   reset();

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[element->index()+offset_] = tmp[element->index()] * element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  DENSESUBVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Subvector operators */
//@{
template< typename VT, bool TF >
inline void reset( DenseSubvector<VT,TF>& dv );

template< typename VT, bool TF >
inline void clear( DenseSubvector<VT,TF>& dv );

template< typename VT, bool TF >
inline bool isDefault( const DenseSubvector<VT,TF>& dv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dense subvector.
// \ingroup dense_subvector
//
// \param dv The dense subvector to be resetted.
// \return void
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void reset( DenseSubvector<VT,TF>& dv )
{
   dv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dense subvector.
// \ingroup dense_subvector
//
// \param dv The dense subvector to be cleared.
// \return void
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void clear( DenseSubvector<VT,TF>& dv )
{
   dv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense subvector is in default state.
// \ingroup dense_subvector
//
// \param dv The dense subvector to be tested for its default state.
// \return \a true in case the given subvector is component-wise zero, \a false otherwise.
//
// This function checks whether the dense subvector is in default state. For instance, in case
// the subvector is instantiated for a vector of built-in integral or floating point data type,
// the function returns \a true in case all subvector elements are 0 and \a false in case any subvector
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isDefault( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline bool isDefault( const DenseSubvector<VT,TF>& dv )
{
   for( size_t i=0UL; i<dv.size(); ++i )
      if( !isDefault( dv[i] ) ) return false;
   return true;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given dense vector.
// \ingroup views
//
// \param dv The dense vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the dense vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given dense
// vector. The following example demonstrates the creation of a subvector of size 8 starting
// from index 4:

   \code
   using blaze::columnVector;

   typedef blaze::DynamicVector<double,columnVector>  Vector;

   Vector v;
   // ... Resizing and initialization
   blaze::DenseSubvector<Vector> = subvector( v, 4UL, 8UL );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DisableIf< Or< IsComputation<VT>, IsTransExpr<VT> >, DenseSubvector<VT> >::Type
   subvector( DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return DenseSubvector<VT>( ~dv, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given dense vector.
// \ingroup views
//
// \param dv The dense vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the dense vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given dense
// vector. The following example demonstrates the creation of a subvector of size 8 starting
// from index 4:

   \code
   using blaze::columnVector;

   typedef blaze::DynamicVector<double,columnVector>  Vector;

   Vector v;
   // ... Resizing and initialization
   blaze::DenseSubvector<Vector> = subvector( v, 4UL, 8UL );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DisableIf< Or< IsComputation<VT>, IsTransExpr<VT> >, DenseSubvector<const VT> >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return DenseSubvector<const VT>( ~dv, index, size );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of another dense subvector.
// \ingroup views
//
// \param dv The constant dense subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the other dense subvector.
//
// This function returns an expression representing the specified subvector of the given
// dense subvector.
*/
template< typename VT  // Type of the dense subvector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,TF>
   subvector( const DenseSubvector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return DenseSubvector<VT,TF>( dv.vector_, dv.offset_ + index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector addition.
// \ingroup views
//
// \param dv The constant vector/vector addition.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the addition.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector addition.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecVecAddExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~dv).leftOperand(), index, size ) + subvector( (~dv).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector subtraction.
// \ingroup views
//
// \param dv The constant vector/vector subtraction.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the subtraction.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector subtraction.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecVecSubExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~dv).leftOperand(), index, size ) - subvector( (~dv).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector multiplication.
// \ingroup views
//
// \param dv The constant vector/vector multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector multiplication.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecVecMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~dv).leftOperand(), index, size ) * subvector( (~dv).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector cross product.
// \ingroup views
//
// \param dv The constant vector/vector cross product.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector cross product.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsCrossExpr<VT>, DynamicVector<typename VT::ElementType,TF> >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   DynamicVector<typename VT::ElementType,TF> tmp( size );

   for( size_t i=0UL; i<size; ++i ) {
      tmp[i] = (~dv)[index+i];
   }

   return tmp;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given matrix/vector multiplication.
// \ingroup views
//
// \param dv The constant matrix/vector multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// matrix/vector multiplication.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsMatVecMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   typename VT::LeftOperand  left ( (~dv).leftOperand()  );
   typename VT::RightOperand right( (~dv).rightOperand() );

   return submatrix( left, index, 0UL, size, left.columns() ) * right;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/matrix multiplication.
// \ingroup views
//
// \param dv The constant vector/matrix multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/matrix multiplication.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsTVecMatMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   typename VT::LeftOperand  left ( (~dv).leftOperand()  );
   typename VT::RightOperand right( (~dv).rightOperand() );

   return left * submatrix( right, 0UL, index, right.rows(), size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar multiplication.
// \ingroup views
//
// \param dv The constant vector/scalar multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar multiplication.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecScalarMultExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~dv).leftOperand(), index, size ) * (~dv).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar division.
// \ingroup views
//
// \param dv The constant vector/scalar division.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar division.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecScalarDivExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~dv).leftOperand(), index, size ) / (~dv).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector abs operation.
// \ingroup views
//
// \param dv The constant vector abs operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the abs operation.
//
// This function returns an expression representing the specified subvector of the given vector
// abs operation.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecAbsExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return abs( subvector( (~dv).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector evaluation operation.
// \ingroup views
//
// \param dv The constant vector evaluation operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the evaluation operation.
//
// This function returns an expression representing the specified subvector of the given vector
// evaluation operation.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecEvalExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return eval( subvector( (~dv).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector transpose operation.
// \ingroup views
//
// \param dv The constant vector transpose operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the transpose operation.
//
// This function returns an expression representing the specified subvector of the given vector
// transpose operation.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsVecTransExpr<VT>, typename SubvectorExprTrait<VT>::Type >::Type
   subvector( const DenseVector<VT,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return trans( subvector( (~dv).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF >
struct SubvectorTrait< DenseSubvector<VT,TF> >
{
   typedef typename SubvectorTrait< typename DenseSubvector<VT,TF>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
