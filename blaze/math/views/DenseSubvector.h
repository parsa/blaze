//=================================================================================================
/*!
//  \file blaze/math/views/DenseSubvector.h
//  \brief Header file for the DenseSubvector class template
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

#ifndef _BLAZE_MATH_VIEWS_DENSESUBVECTOR_H_
#define _BLAZE_MATH_VIEWS_DENSESUBVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Subvector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/CrossExpr.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/Subvector.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Null.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsNumeric.h>
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
/*!\defgroup dense_subvector Dense Subvector
// \ingroup views
*/
/*!\brief View on a specific subvector of a dense vector.
// \ingroup dense_subvector
//
// The DenseSubvector template represents a view on a specific subvector of a dense vector
// primitive. The type of the dense vector is specified via the first template parameter:

   \code
   template< typename VT, bool AF, bool TF >
   class DenseSubvector;
   \endcode

//  - VT: specifies the type of the dense vector primitive. DenseSubvector can be used with every
//        dense vector primitive or view, but does not work with any vector expression type.
//  - AF: the alignment flag specifies whether the subvector is aligned (\a blaze::aligned) or
//        unaligned (\a blaze::unaligned). The default value is \a blaze::unaligned.
//  - TF: specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//        vector (\a blaze::columnVector). This template parameter doesn't have to be explicitly
//        defined, but is automatically derived from the first template parameter.
//
//
// \n \section dense_subvector_setup Setup of Dense Subvectors
//
// A view on a dense subvector can be created very conveniently via the \c subvector() function:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>  DenseVectorType;

   DenseVectorType x;
   // ... Resizing and initialization

   // Create a subvector from index 8 with a size of 16 (i.e. in the range [8..23])
   blaze::DenseSubvector<DenseVectorType> sv = subvector( x, 8UL, 16UL );
   \endcode

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

   // Create a subvector from index 0 with a size of 10 (i.e. in the range [0..9])
   blaze::DenseSubvector<DenseVectorType> sv = subvector( x, 0UL, 10UL );

   // Setting the first ten elements of x to the 2nd row of matrix A
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

// \n \section dense_subvector_arithmetic_operations Arithmetic Operations
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

// \n \section dense_subvector_aligned_subvector Aligned Subvectors
//
// Usually subvectors can be defined anywhere within a vector. They may start at any position and
// may have an arbitrary size (only restricted by the size of the underlying vector). However, in
// contrast to vectors themselves, which are always properly aligned in memory and therefore can
// provide maximum performance, this means that subvectors in general have to be considered to be
// unaligned. This can be made explicit by the \a blaze::unaligned flag:

   \code
   using blaze::unaligned;

   typedef blaze::DynamicVector<double,blaze::rowVector>  DenseVectorType;

   DenseVectorType x;
   // ... Resizing and initialization

   // Identical creations of an unaligned subvector in the range [8..23]
   blaze::DenseSubvector<DenseVectorType>           sv1 = subvector           ( x, 8UL, 16UL );
   blaze::DenseSubvector<DenseVectorType>           sv2 = subvector<unaligned>( x, 8UL, 16UL );
   blaze::DenseSubvector<DenseVectorType,unaligned> sv3 = subvector           ( x, 8UL, 16UL );
   blaze::DenseSubvector<DenseVectorType,unaligned> sv4 = subvector<unaligned>( x, 8UL, 16UL );
   \endcode

// All of these calls to the \c subvector() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned subvector. Whereas this may provide
// full flexibility in the creation of subvectors, this might result in performance restrictions
// (even in case the specified subvector could be aligned). However, it is also possible to create
// aligned subvectors. Aligned subvectors are identical to unaligned subvectors in all aspects,
// except that they may pose additional alignment restrictions and therefore have less flexibility
// during creation, but don't suffer from performance penalties and provide the same performance
// as the underlying vector. Aligned subvectors are created by explicitly specifying the
// \a blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned subvector in the range [8..23]
   blaze::DenseSubvector<DenseVectorType,aligned> sv = subvector<aligned>( x, 8UL, 16UL );
   \endcode

// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). The following source code gives some
// examples for a double precision dense vector, assuming that AVX is available, which packs 4
// \c double values into an intrinsic vector:

   \code
   using blaze::columnVector;

   typedef blaze::DynamicVector<double,columnVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType,aligned>  SubvectorType;

   VectorType d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning and the size is a multiple of 4
   SubvectorType dsv1 = subvector<aligned>( d, 0UL, 12UL );

   // OK: Start index and the size are both a multiple of 4
   SubvectorType dsv2 = subvector<aligned>( d, 4UL, 8UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   SubvectorType dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4
   SubvectorType dsv4 = subvector<aligned>( d, 5UL, 8UL );

   // Error: Size is not a multiple of 4 and the subvector does not include the last element
   SubvectorType dsv5 = subvector<aligned>( d, 8UL, 5UL );
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
template< typename VT                         // Type of the dense vector
        , bool AF = unaligned                 // Alignment flag
        , bool TF = IsRowVector<VT>::value >  // Transpose flag
class DenseSubvector : public DenseVector< DenseSubvector<VT,AF,TF>, TF >
                     , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense vector expression.
   typedef typename If< IsExpression<VT>, VT, VT& >::Type  Operand;

   //! Intrinsic trait for the vector element type.
   typedef IntrinsicTrait<typename VT::ElementType>  IT;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<VT,AF,TF>            This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<VT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename VT::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the subvector elements.
   typedef typename VT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubvector&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   typedef typename VT::ConstReference  ConstReference;

   //! Reference to a non-constant subvector value.
   typedef typename If< IsConst<VT>, ConstReference, typename VT::Reference >::Type  Reference;

   //! Pointer to a constant subvector value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant subvector value.
   typedef typename If< Or< IsConst<VT>, Not< HasMutableDataAccess<VT> > >
                      , ConstPointer, ElementType* >::Type  Pointer;
   //**********************************************************************************************

   //**SubvectorIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse subvector.
   */
   template< typename IteratorType >  // Type of the dense vector iterator
   class SubvectorIterator
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
      /*!\brief Default constructor of the SubvectorIterator class.
      */
      inline SubvectorIterator()
         : iterator_ (       )  // Iterator to the current subvector element
         , isAligned_( false )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the SubvectorIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param isMemoryAligned Memory alignment flag.
      */
      inline SubvectorIterator( IteratorType iterator, bool isMemoryAligned )
         : iterator_ ( iterator        )  // Iterator to the current subvector element
         , isAligned_( isMemoryAligned )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubvectorIterator instances.
      //
      // \param it The subvector iterator to be copied
      */
      template< typename IteratorType2 >
      inline SubvectorIterator( const SubvectorIterator<IteratorType2>& it )
         : iterator_ ( it.base()      )  // Iterator to the current subvector element
         , isAligned_( it.isAligned() )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline SubvectorIterator& operator+=( size_t inc ) {
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
      inline SubvectorIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubvectorIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubvectorIterator operator++( int ) {
         return SubvectorIterator( iterator_++, isAligned_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline SubvectorIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubvectorIterator operator--( int ) {
         return SubvectorIterator( iterator_--, isAligned_ );
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
      /*!\brief Aligned load of an intrinsic element of the dense subvector.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an aligned load of the current intrinsic element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType load() const {
         return loadu();
      }
      //*******************************************************************************************

      //**Loada function***************************************************************************
      /*!\brief Aligned load of an intrinsic element of the dense subvector.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an aligned load of the current intrinsic element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType loada() const {
         return iterator_.loada();
      }
      //*******************************************************************************************

      //**Loadu function***************************************************************************
      /*!\brief Unaligned load of an intrinsic element of the dense subvector.
      //
      // \return The loaded intrinsic element.
      //
      // This function performs an unaligned load of the current intrinsic element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline IntrinsicType loadu() const {
         if( isAligned_ ) {
            return iterator_.loada();
         }
         else {
            return iterator_.loadu();
         }
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const SubvectorIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const SubvectorIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const SubvectorIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const SubvectorIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const SubvectorIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const SubvectorIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const SubvectorIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a SubvectorIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const SubvectorIterator operator+( const SubvectorIterator& it, size_t inc ) {
         return SubvectorIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a SubvectorIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const SubvectorIterator operator+( size_t inc, const SubvectorIterator& it ) {
         return SubvectorIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a SubvectorIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const SubvectorIterator operator-( const SubvectorIterator& it, size_t dec ) {
         return SubvectorIterator( it.iterator_ - dec, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the subvector iterator.
      //
      // \return The current position of the subvector iterator.
      */
      inline IteratorType base() const {
         return iterator_;
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
      IteratorType iterator_;   //!< Iterator to the current subvector element.
      bool         isAligned_;  //!< Memory alignment flag.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubvectorIterator<typename VT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<VT>, ConstIterator, SubvectorIterator<typename VT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT::vectorizable };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = VT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubvector( Operand vector, size_t index, size_t n );
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
                            inline DenseSubvector& operator= ( const ElementType& rhs );
                            inline DenseSubvector& operator= ( const DenseSubvector& rhs );
   template< typename VT2 > inline DenseSubvector& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator*=( const DenseVector<VT2,TF>&  rhs );
   template< typename VT2 > inline DenseSubvector& operator*=( const SparseVector<VT2,TF>& rhs );

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
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedAddAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
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
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
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
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
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

   template< typename VT2, bool AF2, bool TF2 >
   inline bool canAlias( const DenseSubvector<VT2,AF2,TF2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename VT2, bool AF2, bool TF2 >
   inline bool isAliased( const DenseSubvector<VT2,AF2,TF2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load ( size_t index ) const;
   inline IntrinsicType loada( size_t index ) const;
   inline IntrinsicType loadu( size_t index ) const;

   inline void store ( size_t index, const IntrinsicType& value );
   inline void storea( size_t index, const IntrinsicType& value );
   inline void storeu( size_t index, const IntrinsicType& value );
   inline void stream( size_t index, const IntrinsicType& value );

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
   const bool isAligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the subvector is fully aligned
                                with respect to the given element type and the available instruction
                                set. In case the subvector is fully aligned it is possible to use
                                aligned loads and stores instead of unaligned loads and stores. In
                                order to be aligned, the first element of the subvector must be
                                aligned. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename VT2, bool AF2, bool TF2 > friend class DenseSubvector;

   template< bool AF1, typename VT2, bool AF2, bool TF2 >
   friend const DenseSubvector<VT2,AF1,TF2>
      subvector( const DenseSubvector<VT2,AF2,TF2>& dv, size_t index, size_t size );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isIntact( const DenseSubvector<VT2,AF2,TF2>& dv );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const DenseSubvector<VT2,AF2,TF2>& a, const DenseVector<VT2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const DenseVector<VT2,TF2>& a, const DenseSubvector<VT2,AF2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const DenseSubvector<VT2,AF2,TF2>& a, const DenseSubvector<VT2,AF2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryAddAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool trySubAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryMultAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2 >
   friend typename DerestrictTrait< DenseSubvector<VT2,AF2,TF2> >::Type
      derestrict( DenseSubvector<VT2,AF2,TF2>& dv );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBVECTOR_TYPE  ( VT );
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
//
// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the size of the given vector or the subvector is specified beyond the size of the vector)
// a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,AF,TF>::DenseSubvector( Operand vector, size_t index, size_t n )
   : vector_   ( vector )  // The vector containing the subvector
   , offset_   ( index  )  // The offset of the subvector within the dense vector
   , size_     ( n      )  // The size of the subvector
   , isAligned_( vectorizable && vector.data() != NULL && checkAlignment( data() ) )
{
   if( index + n > vector.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }
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
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::Reference
   DenseSubvector<VT,AF,TF>::operator[]( size_t index )
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
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstReference
   DenseSubvector<VT,AF,TF>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return const_cast<const VT&>( vector_ )[offset_+index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::Reference
   DenseSubvector<VT,AF,TF>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstReference
   DenseSubvector<VT,AF,TF>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::Pointer DenseSubvector<VT,AF,TF>::data()
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstPointer DenseSubvector<VT,AF,TF>::data() const
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::Iterator DenseSubvector<VT,AF,TF>::begin()
{
   return Iterator( vector_.begin() + offset_, isAligned_ );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstIterator DenseSubvector<VT,AF,TF>::begin() const
{
   return ConstIterator( vector_.cbegin() + offset_, isAligned_ );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstIterator DenseSubvector<VT,AF,TF>::cbegin() const
{
   return ConstIterator( vector_.cbegin() + offset_, isAligned_ );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::Iterator DenseSubvector<VT,AF,TF>::end()
{
   return Iterator( vector_.begin() + offset_ + size_, isAligned_ );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstIterator DenseSubvector<VT,AF,TF>::end() const
{
   return ConstIterator( vector_.cbegin() + offset_ + size_, isAligned_ );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::ConstIterator DenseSubvector<VT,AF,TF>::cend() const
{
   return ConstIterator( vector_.cbegin() + offset_ + size_, isAligned_ );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,AF,TF>& DenseSubvector<VT,AF,TF>::operator=( const ElementType& rhs )
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
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,AF,TF>& DenseSubvector<VT,AF,TF>::operator=( const DenseSubvector& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( &rhs == this || ( &vector_ == &rhs.vector_ && offset_ == rhs.offset_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Subvector sizes do not match" );
   }

   if( !tryAssign( vector_, rhs, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &vector_ ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,AF,TF>& DenseSubvector<VT,AF,TF>::operator=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector<VT2>::value )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,AF,TF>& DenseSubvector<VT,AF,TF>::operator+=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !tryAddAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,AF,TF>& DenseSubvector<VT,AF,TF>::operator-=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !trySubAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a dense vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline DenseSubvector<VT,AF,TF>&
   DenseSubvector<VT,AF,TF>::operator*=( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !tryMultAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side sparse vector
inline DenseSubvector<VT,AF,TF>&
   DenseSubvector<VT,AF,TF>::operator*=( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const ResultType tmp( *this * (~rhs) );

   if( !tryAssign( vector_, tmp, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

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
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubvector<VT,AF,TF> >::Type&
   DenseSubvector<VT,AF,TF>::operator*=( Other rhs )
{
   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) * rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a subvector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the assigned subvector.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename VT       // Type of the dense vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubvector<VT,AF,TF> >::Type&
   DenseSubvector<VT,AF,TF>::operator/=( Other rhs )
{
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
/*!\brief Returns the current size/dimension of the dense subvector.
//
// \return The size of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,AF,TF>::size() const
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,AF,TF>::capacity() const
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,AF,TF>::nonZeros() const
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,AF,TF>::reset()
{
   using blaze::clear;

   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i )
      clear( vector_[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the dense subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the scalar value
inline DenseSubvector<VT,AF,TF>& DenseSubvector<VT,AF,TF>::scale( const Other& scalar )
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
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubvector<VT,AF,TF>::canAlias( const Other* alias ) const
{
   return vector_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense subvector can alias with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT   // Type of the dense vector
        , bool AF       // Alignment flag
        , bool TF >     // Transpose flag
template< typename VT2  // Data type of the foreign dense subvector
        , bool AF2      // Alignment flag of the foreign dense subvector
        , bool TF2 >    // Transpose flag of the foreign dense subvector
inline bool DenseSubvector<VT,AF,TF>::canAlias( const DenseSubvector<VT2,AF2,TF2>* alias ) const
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset_ + size_ > alias->offset_ ) && ( offset_ < alias->offset_ + alias->size_ ) );
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
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool DenseSubvector<VT,AF,TF>::isAliased( const Other* alias ) const
{
   return vector_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the dense subvector is aliased with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT   // Type of the dense vector
        , bool AF       // Alignment flag
        , bool TF >     // Transpose flag
template< typename VT2  // Data type of the foreign dense subvector
        , bool AF2      // Alignment flag of the foreign dense subvector
        , bool TF2 >    // Transpose flag of the foreign dense subvector
inline bool DenseSubvector<VT,AF,TF>::isAliased( const DenseSubvector<VT2,AF2,TF2>* alias ) const
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset_ + size_ > alias->offset_ ) && ( offset_ < alias->offset_ + alias->size_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the subvector is properly aligned in memory.
//
// \return \a true in case the subvector is aligned, \a false if not.
//
// This function returns whether the subvector is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the subvector are guaranteed to conform to the
// alignment restrictions of the underlying element type.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool DenseSubvector<VT,AF,TF>::isAligned() const
{
   return isAligned_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the subvector can be used in SMP assignments.
//
// \return \a true in case the subvector can be used in SMP assignments, \a false if not.
//
// This function returns whether the subvector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// subvector).
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool DenseSubvector<VT,AF,TF>::canSMPAssign() const
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
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
template< typename VT       // Type of the dense vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::IntrinsicType
   DenseSubvector<VT,AF,TF>::load( size_t index ) const
{
   if( isAligned_ )
      return loada( index );
   else
      return loadu( index );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::IntrinsicType
   DenseSubvector<VT,AF,TF>::loada( size_t index ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   return vector_.loada( offset_+index );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,AF,TF>::IntrinsicType
   DenseSubvector<VT,AF,TF>::loadu( size_t index ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   return vector_.loadu( offset_+index );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,AF,TF>::store( size_t index, const IntrinsicType& value )
{
   if( isAligned_ )
      storea( index, value );
   else
      storeu( index, value );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,AF,TF>::storea( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   vector_.storea( offset_+index, value );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,AF,TF>::storeu( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   vector_.storeu( offset_+index, value );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,AF,TF>::stream( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   if( isAligned_ )
      vector_.stream( offset_+index, value );
   else
      vector_.storeu( offset_+index, value );
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] = (~rhs)[i    ];
      vector_[offset_+i+1UL] = (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] = (~rhs)[ipos];
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   if( useStreaming && isAligned_ &&
       ( size_ > ( cacheSize/( sizeof(ElementType) * 3UL ) ) ) &&
       !(~rhs).isAliased( &vector_ ) )
   {
      size_t i( 0UL );

      for( ; i<ipos; i+=IT::size ) {
         vector_.stream( offset_+i, (~rhs).load(i) );
      }
      for( ; i<size_; ++i ) {
         vector_[offset_+i] = (~rhs)[i];
      }
   }
   else
   {
      size_t i( 0UL );
      typename VT2::ConstIterator it( (~rhs).begin() );

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
         store( i             , it.load() ); it += IT::size;
         store( i+IT::size    , it.load() ); it += IT::size;
         store( i+IT::size*2UL, it.load() ); it += IT::size;
         store( i+IT::size*3UL, it.load() ); it += IT::size;
      }
      for( ; i<ipos; i+=IT::size, it+=IT::size ) {
         store( i, it.load() );
      }
      for( ; i<size_; ++i, ++it ) {
         vector_[offset_+i] = *it;
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,AF,TF>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] = element->value();
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedAddAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] += (~rhs)[i    ];
      vector_[offset_+i+1UL] += (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] += (~rhs)[ipos];
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedAddAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT2::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      store( i             , load(i             ) + it.load() ); it += IT::size;
      store( i+IT::size    , load(i+IT::size    ) + it.load() ); it += IT::size;
      store( i+IT::size*2UL, load(i+IT::size*2UL) + it.load() ); it += IT::size;
      store( i+IT::size*3UL, load(i+IT::size*3UL) + it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      store( i, load(i) + it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      vector_[offset_+i] += *it;
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,AF,TF>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] += element->value();
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedSubAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] -= (~rhs)[i    ];
      vector_[offset_+i+1UL] -= (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] -= (~rhs)[ipos];
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedSubAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT2::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      store( i             , load(i             ) - it.load() ); it += IT::size;
      store( i+IT::size    , load(i+IT::size    ) - it.load() ); it += IT::size;
      store( i+IT::size*2UL, load(i+IT::size*2UL) - it.load() ); it += IT::size;
      store( i+IT::size*3UL, load(i+IT::size*3UL) - it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      store( i, load(i) - it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      vector_[offset_+i] -= *it;
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,AF,TF>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] -= element->value();
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedMultAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] *= (~rhs)[i    ];
      vector_[offset_+i+1UL] *= (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] *= (~rhs)[ipos];
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,AF,TF>::BLAZE_TEMPLATE VectorizedMultAssign<VT2> >::Type
   DenseSubvector<VT,AF,TF>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT2::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      store( i             , load(i             ) * it.load() ); it += IT::size;
      store( i+IT::size    , load(i+IT::size    ) * it.load() ); it += IT::size;
      store( i+IT::size*2UL, load(i+IT::size*2UL) * it.load() ); it += IT::size;
      store( i+IT::size*3UL, load(i+IT::size*3UL) * it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      store( i, load(i) * it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      vector_[offset_+i] *= *it;
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
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,AF,TF>::multAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] = tmp[element->index()] * element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ALIGNED SUBVECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubvector for aligned subvectors.
// \ingroup dense_subvector
//
// This specialization of DenseSubvector adapts the class template to the requirements of
// aligned subvectors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DenseSubvector<VT,aligned,TF> : public DenseVector< DenseSubvector<VT,aligned,TF>, TF >
                                    , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense vector expression.
   typedef typename If< IsExpression<VT>, VT, VT& >::Type  Operand;

   //! Intrinsic trait for the vector element type.
   typedef IntrinsicTrait<typename VT::ElementType>  IT;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<VT,aligned,TF>       This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<VT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename VT::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename IT::Type                   IntrinsicType;  //!< Intrinsic type of the subvector elements.
   typedef typename VT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const DenseSubvector&               CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   typedef typename VT::ConstReference  ConstReference;

   //! Reference to a non-constant subvector value.
   typedef typename If< IsConst<VT>, ConstReference, typename VT::Reference >::Type  Reference;

   //! Pointer to a constant subvector value.
   typedef const ElementType*  ConstPointer;

   //! Pointer to a non-constant subvector value.
   typedef typename If< Or< IsConst<VT>, Not< HasMutableDataAccess<VT> > >
                      , ConstPointer, ElementType* >::Type  Pointer;

   //! Iterator over constant elements.
   typedef typename VT::ConstIterator  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<VT>, ConstIterator, typename VT::Iterator >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT::vectorizable };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = VT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DenseSubvector( Operand vector, size_t index, size_t n );
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
                            inline DenseSubvector& operator= ( const ElementType& rhs );
                            inline DenseSubvector& operator= ( const DenseSubvector& rhs );
   template< typename VT2 > inline DenseSubvector& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline DenseSubvector& operator*=( const DenseVector<VT2,TF>&  rhs );
   template< typename VT2 > inline DenseSubvector& operator*=( const SparseVector<VT2,TF>& rhs );

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
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedAddAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::addition };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedSubAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::subtraction };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct VectorizedMultAssign {
      enum { value = useOptimizedKernels &&
                     vectorizable && VT2::vectorizable &&
                     IsSame<ElementType,typename VT2::ElementType>::value &&
                     IntrinsicTrait<ElementType>::multiplication };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const;

   template< typename VT2, bool AF2, bool TF2 >
   inline bool canAlias( const DenseSubvector<VT2,AF2,TF2>* alias ) const;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const;

   template< typename VT2, bool AF2, bool TF2 >
   inline bool isAliased( const DenseSubvector<VT2,AF2,TF2>* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load ( size_t index ) const;
   inline IntrinsicType loada( size_t index ) const;
   inline IntrinsicType loadu( size_t index ) const;

   inline void store ( size_t index, const IntrinsicType& value );
   inline void storea( size_t index, const IntrinsicType& value );
   inline void storeu( size_t index, const IntrinsicType& value );
   inline void stream( size_t index, const IntrinsicType& value );

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
   Operand      vector_;  //!< The dense vector containing the subvector.
   const size_t offset_;  //!< The offset of the subvector within the dense vector.
   const size_t size_;    //!< The size of the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename VT2, bool AF2, bool TF2 > friend class DenseSubvector;

   template< bool AF1, typename VT2, bool AF2, bool TF2 >
   friend const DenseSubvector<VT2,AF1,TF2>
      subvector( const DenseSubvector<VT2,AF2,TF2>& dv, size_t index, size_t size );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isIntact( const DenseSubvector<VT2,AF2,TF2>& dv );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const DenseSubvector<VT2,AF2,TF2>& a, const DenseVector<VT2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const DenseVector<VT2,TF2>& a, const DenseSubvector<VT2,AF2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const DenseSubvector<VT2,AF2,TF2>& a, const DenseSubvector<VT2,AF2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryAddAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool trySubAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryMultAssign( const DenseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2 >
   friend typename DerestrictTrait< DenseSubvector<VT2,AF2,TF2> >::Type
      derestrict( DenseSubvector<VT2,AF2,TF2>& dv );
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBVECTOR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
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
/*!\brief The constructor for DenseSubvector.
//
// \param vector The dense vector containing the subvector.
// \param index The first index of the subvector in the given vector.
// \param n The size of the subvector.
// \exception std::invalid_argument Invalid subvector specification.
//
// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the size of the given vector or the subvector is specified beyond the size of the vector)
// a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,aligned,TF>::DenseSubvector( Operand vector, size_t index, size_t n )
   : vector_( vector )  // The vector containing the subvector
   , offset_( index  )  // The offset of the subvector within the dense vector
   , size_  ( n      )  // The size of the subvector
{
   if( index + n > vector.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }

   if( vectorizable && vector_.data() != NULL && !checkAlignment( data() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector alignment" );
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
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::Reference
   DenseSubvector<VT,aligned,TF>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return vector_[offset_+index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstReference
   DenseSubvector<VT,aligned,TF>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return const_cast<const VT&>( vector_ )[offset_+index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::Reference
   DenseSubvector<VT,aligned,TF>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstReference
   DenseSubvector<VT,aligned,TF>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::Pointer DenseSubvector<VT,aligned,TF>::data()
{
   return vector_.data() + offset_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstPointer DenseSubvector<VT,aligned,TF>::data() const
{
   return vector_.data() + offset_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::Iterator DenseSubvector<VT,aligned,TF>::begin()
{
   return ( vector_.begin() + offset_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstIterator
   DenseSubvector<VT,aligned,TF>::begin() const
{
   return ( vector_.cbegin() + offset_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstIterator
   DenseSubvector<VT,aligned,TF>::cbegin() const
{
   return ( vector_.cbegin() + offset_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::Iterator DenseSubvector<VT,aligned,TF>::end()
{
   return ( vector_.begin() + offset_ + size_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstIterator
   DenseSubvector<VT,aligned,TF>::end() const
{
   return ( vector_.cbegin() + offset_ + size_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline typename DenseSubvector<VT,aligned,TF>::ConstIterator
   DenseSubvector<VT,aligned,TF>::cend() const
{
   return ( vector_.cbegin() + offset_ + size_ );
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
/*!\brief Homogenous assignment to all subvector elements.
//
// \param rhs Scalar value to be assigned to all subvector elements.
// \return Reference to the assigned subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator=( const ElementType& rhs )
{
   const size_t iend( offset_ + size_ );

   for( size_t i=offset_; i<iend; ++i )
      vector_[i] = rhs;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DenseSubvector.
//
// \param rhs Dense subvector to be copied.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Subvector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator=( const DenseSubvector& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( &rhs == this || ( &vector_ == &rhs.vector_ && offset_ == rhs.offset_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Subvector sizes do not match" );
   }

   if( !tryAssign( vector_, rhs, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &vector_ ) ) {
      const ResultType tmp( ~rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !tryAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector<VT2>::value )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator+=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !tryAddAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator-=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !trySubAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a dense vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator*=( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( typename VT2::ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   typedef typename If< IsRestricted<VT>, typename VT2::CompositeType, const VT2& >::Type  Right;
   Right right( ~rhs );

   if( !tryMultAssign( vector_, right, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( IsReference<Right>::value && right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side sparse vector
inline DenseSubvector<VT,aligned,TF>&
   DenseSubvector<VT,aligned,TF>::operator*=( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const ResultType tmp( *this * (~rhs) );

   if( !tryAssign( vector_, tmp, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a subvector and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the assigned subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubvector<VT,aligned,TF> >::Type&
   DenseSubvector<VT,aligned,TF>::operator*=( Other rhs )
{
   typename DerestrictTrait<This>::Type left( derestrict( *this ) );
   smpAssign( left, (*this) * rhs );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a subvector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the assigned subvector.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DenseSubvector<VT,aligned,TF> >::Type&
   DenseSubvector<VT,aligned,TF>::operator/=( Other rhs )
{
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
/*!\brief Returns the current size/dimension of the dense subvector.
//
// \return The size of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,aligned,TF>::size() const
{
   return size_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense subvector.
//
// \return The capacity of the dense subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,aligned,TF>::capacity() const
{
   return vector_.capacity() - offset_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the subvector.
//
// \return The number of non-zero elements in the subvector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the subvector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline size_t DenseSubvector<VT,aligned,TF>::nonZeros() const
{
   size_t nonzeros( 0 );

   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i ) {
      if( !isDefault( vector_[i] ) )
         ++nonzeros;
   }

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
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void DenseSubvector<VT,aligned,TF>::reset()
{
   using blaze::clear;

   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i )
      clear( vector_[i] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the scalar value
inline DenseSubvector<VT,aligned,TF>& DenseSubvector<VT,aligned,TF>::scale( const Other& scalar )
{
   const size_t iend( offset_ + size_ );
   for( size_t i=offset_; i<iend; ++i )
      vector_[i] *= scalar;
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
inline bool DenseSubvector<VT,aligned,TF>::canAlias( const Other* alias ) const
{
   return vector_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector can alias with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT   // Type of the dense vector
        , bool TF >     // Transpose flag
template< typename VT2  // Data type of the foreign dense subvector
        , bool AF2      // Alignment flag of the foreign dense subvector
        , bool TF2 >    // Transpose flag of the foreign dense subvector
inline bool DenseSubvector<VT,aligned,TF>::canAlias( const DenseSubvector<VT2,AF2,TF2>* alias ) const
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset_ + size_ > alias->offset_ ) && ( offset_ < alias->offset_ + alias->size_ ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
inline bool DenseSubvector<VT,aligned,TF>::isAliased( const Other* alias ) const
{
   return vector_.isAliased( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector is aliased with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT   // Type of the dense vector
        , bool TF >     // Transpose flag
template< typename VT2  // Data type of the foreign dense subvector
        , bool AF2      // Alignment flag of the foreign dense subvector
        , bool TF2 >    // Transpose flag of the foreign dense subvector
inline bool DenseSubvector<VT,aligned,TF>::isAliased( const DenseSubvector<VT2,AF2,TF2>* alias ) const
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset_ + size_ > alias->offset_ ) && ( offset_ < alias->offset_ + alias->size_ ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the subvector is properly aligned in memory.
//
// \return \a true in case the subvector is aligned, \a false if not.
//
// This function returns whether the subvector is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the subvector are guaranteed to conform to the
// alignment restrictions of the underlying element type.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline bool DenseSubvector<VT,aligned,TF>::isAligned() const
{
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the subvector can be used in SMP assignments.
//
// \return \a true in case the subvector can be used in SMP assignments, \a false if not.
//
// This function returns whether the subvector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// subvector).
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline bool DenseSubvector<VT,aligned,TF>::canSMPAssign() const
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of an intrinsic element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded intrinsic element.
//
// This function performs a load of a specific intrinsic element of the dense subvector. The
// index must be smaller than the number of subvector elements and it must be a multiple of
// the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
BLAZE_ALWAYS_INLINE typename DenseSubvector<VT,aligned,TF>::IntrinsicType
   DenseSubvector<VT,aligned,TF>::load( size_t index ) const
{
   return loada( index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename VT       // Type of the dense vector
        , bool TF >         // Transpose flag
BLAZE_ALWAYS_INLINE typename DenseSubvector<VT,aligned,TF>::IntrinsicType
   DenseSubvector<VT,aligned,TF>::loada( size_t index ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   return vector_.loada( offset_+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
BLAZE_ALWAYS_INLINE typename DenseSubvector<VT,aligned,TF>::IntrinsicType
   DenseSubvector<VT,aligned,TF>::loadu( size_t index ) const
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   return vector_.loadu( offset_+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of an intrinsic element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs a store a specific intrinsic element of the dense subvector. The
// index must be smaller than the number of subvector elements and it must be a multiple of
// the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
BLAZE_ALWAYS_INLINE void
   DenseSubvector<VT,aligned,TF>::store( size_t index, const IntrinsicType& value )
{
   storea( index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
BLAZE_ALWAYS_INLINE void
   DenseSubvector<VT,aligned,TF>::storea( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   vector_.storea( offset_+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
BLAZE_ALWAYS_INLINE void
   DenseSubvector<VT,aligned,TF>::storeu( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   vector_.storeu( offset_+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
BLAZE_ALWAYS_INLINE void
   DenseSubvector<VT,aligned,TF>::stream( size_t index, const IntrinsicType& value )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + IT::size <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL   , "Invalid subvector access index" );

   vector_.stream( offset_+index, value );
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] = (~rhs)[i    ];
      vector_[offset_+i+1UL] = (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] = (~rhs)[ipos];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   if( useStreaming && size_ > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(~rhs).isAliased( &vector_ ) )
   {
      size_t i( 0UL );

      for( ; i<ipos; i+=IT::size ) {
         stream( i, (~rhs).load(i) );
      }
      for( ; i<size_; ++i ) {
         vector_[offset_+i] = (~rhs)[i];
      }
   }
   else
   {
      size_t i( 0UL );
      typename VT2::ConstIterator it( (~rhs).begin() );

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
         store( i             , it.load() ); it += IT::size;
         store( i+IT::size    , it.load() ); it += IT::size;
         store( i+IT::size*2UL, it.load() ); it += IT::size;
         store( i+IT::size*3UL, it.load() ); it += IT::size;
      }
      for( ; i<ipos; i+=IT::size, it+=IT::size ) {
         store( i, it.load() );
      }
      for( ; i<size_; ++i, ++it ) {
         vector_[offset_+i] = *it;
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,aligned,TF>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] = element->value();
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedAddAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] += (~rhs)[i    ];
      vector_[offset_+i+1UL] += (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] += (~rhs)[ipos];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedAddAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT2::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      store( i             , load(i             ) + it.load() ); it += IT::size;
      store( i+IT::size    , load(i+IT::size    ) + it.load() ); it += IT::size;
      store( i+IT::size*2UL, load(i+IT::size*2UL) + it.load() ); it += IT::size;
      store( i+IT::size*3UL, load(i+IT::size*3UL) + it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      store( i, load(i) + it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      vector_[offset_+i] += *it;
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,aligned,TF>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] += element->value();
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedSubAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] -= (~rhs)[i    ];
      vector_[offset_+i+1UL] -= (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] -= (~rhs)[ipos];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedSubAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT2::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      store( i             , load(i             ) - it.load() ); it += IT::size;
      store( i+IT::size    , load(i+IT::size    ) - it.load() ); it += IT::size;
      store( i+IT::size*2UL, load(i+IT::size*2UL) - it.load() ); it += IT::size;
      store( i+IT::size*3UL, load(i+IT::size*3UL) - it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      store( i, load(i) - it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      vector_[offset_+i] -= *it;
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,aligned,TF>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] -= element->value();
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename DisableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedMultAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size() & size_t(-2) );
   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset_+i    ] *= (~rhs)[i    ];
      vector_[offset_+i+1UL] *= (~rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset_+ipos] *= (~rhs)[ipos];
   }
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline typename EnableIf< typename DenseSubvector<VT,aligned,TF>::BLAZE_TEMPLATE VectorizedMultAssign<VT2> >::Type
   DenseSubvector<VT,aligned,TF>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( size_ & size_t(-IT::size) );
   BLAZE_INTERNAL_ASSERT( ( size_ - ( size_ % (IT::size) ) ) == ipos, "Invalid end calculation" );

   size_t i( 0UL );
   typename VT2::ConstIterator it( (~rhs).begin() );

   for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL ) {
      store( i             , load(i             ) * it.load() ); it += IT::size;
      store( i+IT::size    , load(i+IT::size    ) * it.load() ); it += IT::size;
      store( i+IT::size*2UL, load(i+IT::size*2UL) * it.load() ); it += IT::size;
      store( i+IT::size*3UL, load(i+IT::size*3UL) * it.load() ); it += IT::size;
   }
   for( ; i<ipos; i+=IT::size, it+=IT::size ) {
      store( i, load(i) * it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      vector_[offset_+i] *= *it;
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
template< typename VT     // Type of the dense vector
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void DenseSubvector<VT,aligned,TF>::multAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      vector_[offset_+element->index()] = tmp[element->index()] * element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DVECDVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubvector for dense vector/dense vector cross products.
// \ingroup dense_subvector
//
// This specialization of DenseSubvector adapts the class template to the special case of
// dense vector/dense vector cross products.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
class DenseSubvector< DVecDVecCrossExpr<VT1,VT2>, unaligned, false >
   : public DenseVector< DenseSubvector< DVecDVecCrossExpr<VT1,VT2>, unaligned, false >, false >
   , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   typedef DVecDVecCrossExpr<VT1,VT2>  CPE;  //!< Type of the cross product expression.
   typedef typename CPE::ResultType    RT;   //!< Result type of the cross product expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<CPE,unaligned,false>  This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<RT>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType   TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename CPE::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename CPE::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const ResultType                     CompositeType;  //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DenseSubvector specialization class.
   //
   // \param vector The dense vector/dense vector cross product expression.
   // \param index The first index of the subvector in the given expression.
   // \param n The size of the subvector.
   */
   explicit inline DenseSubvector( const CPE& vector, size_t index, size_t n )
      : vector_( vector )  // The dense vector/dense vector cross product expression
      , offset_( index  )  // The offset of the subvector within the cross product expression
      , size_  ( n      )  // The size of the subvector
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset_+index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return size_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const {
      return vector_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CPE          vector_;  //!< The dense vector/dense vector cross product expression.
   const size_t offset_;  //!< The offset of the subvector within the cross product expression.
   const size_t size_;    //!< The size of the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< bool AF1, typename VT, bool AF2, bool TF >
   friend const DenseSubvector<VT,AF1,TF>
      subvector( const DenseSubvector<VT,AF2,TF>& dv, size_t index, size_t size );

   template< typename VT3, bool AF, bool TF >
   friend bool isIntact( const DenseSubvector<VT3,AF,TF>& dv );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseVector<VT3,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseVector<VT3,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAddAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool trySubAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryMultAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DVECSVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubvector for dense vector/sparse vector cross products.
// \ingroup dense_subvector
//
// This specialization of DenseSubvector adapts the class template to the special case of
// dense vector/sparse vector cross products.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side sparse vector
class DenseSubvector< DVecSVecCrossExpr<VT1,VT2>, unaligned, false >
   : public DenseVector< DenseSubvector< DVecSVecCrossExpr<VT1,VT2>, unaligned, false >, false >
   , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   typedef DVecSVecCrossExpr<VT1,VT2>  CPE;  //!< Type of the cross product expression.
   typedef typename CPE::ResultType    RT;   //!< Result type of the cross product expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<CPE,unaligned,false>  This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<RT>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType   TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename CPE::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename CPE::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const ResultType                     CompositeType;  //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DenseSubvector specialization class.
   //
   // \param vector The dense vector/sparse vector cross product expression.
   // \param index The first index of the subvector in the given expression.
   // \param n The size of the subvector.
   */
   explicit inline DenseSubvector( const CPE& vector, size_t index, size_t n )
      : vector_( vector )  // The dense vector/sparse vector cross product expression
      , offset_( index  )  // The offset of the subvector within the cross product expression
      , size_  ( n      )  // The size of the subvector
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset_+index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return size_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const {
      return vector_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CPE          vector_;  //!< The dense vector/sparse vector cross product expression.
   const size_t offset_;  //!< The offset of the subvector within the cross product expression.
   const size_t size_;    //!< The size of the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< bool AF1, typename VT, bool AF2, bool TF >
   friend const DenseSubvector<VT,AF1,TF>
      subvector( const DenseSubvector<VT,AF2,TF>& dv, size_t index, size_t size );

   template< typename VT3, bool AF, bool TF >
   friend bool isIntact( const DenseSubvector<VT3,AF,TF>& dv );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseVector<VT3,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseVector<VT3,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAddAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool trySubAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryMultAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SVECDVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubvector for sparse vector/dense vector cross products.
// \ingroup dense_subvector
//
// This specialization of DenseSubvector adapts the class template to the special case of
// sparse vector/dense vector cross products.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
class DenseSubvector< SVecDVecCrossExpr<VT1,VT2>, unaligned, false >
   : public DenseVector< DenseSubvector< SVecDVecCrossExpr<VT1,VT2>, unaligned, false >, false >
   , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   typedef SVecDVecCrossExpr<VT1,VT2>  CPE;  //!< Type of the cross product expression.
   typedef typename CPE::ResultType    RT;   //!< Result type of the cross product expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<CPE,unaligned,false>  This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<RT>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType   TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename CPE::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename CPE::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const ResultType                     CompositeType;  //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DenseSubvector specialization class.
   //
   // \param vector The sparse vector/dense vector cross product expression.
   // \param index The first index of the subvector in the given expression.
   // \param n The size of the subvector.
   */
   explicit inline DenseSubvector( const CPE& vector, size_t index, size_t n )
      : vector_( vector )  // The sparse vector/dense vector cross product expression
      , offset_( index  )  // The offset of the subvector within the cross product expression
      , size_  ( n      )  // The size of the subvector
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset_+index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return size_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const {
      return vector_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CPE          vector_;  //!< The sparse vector/dense vector cross product expression.
   const size_t offset_;  //!< The offset of the subvector within the cross product expression.
   const size_t size_;    //!< The size of the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< bool AF1, typename VT, bool AF2, bool TF >
   friend const DenseSubvector<VT,AF1,TF>
      subvector( const DenseSubvector<VT,AF2,TF>& dv, size_t index, size_t size );

   template< typename VT3, bool AF, bool TF >
   friend bool isIntact( const DenseSubvector<VT3,AF,TF>& dv );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseVector<VT3,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseVector<VT3,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAddAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool trySubAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryMultAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SVECSVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DenseSubvector for sparse vector/sparse vector cross products.
// \ingroup dense_subvector
//
// This specialization of DenseSubvector adapts the class template to the special case of
// sparse vector/sparse vector cross products.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
class DenseSubvector< SVecSVecCrossExpr<VT1,VT2>, unaligned, false >
   : public DenseVector< DenseSubvector< SVecSVecCrossExpr<VT1,VT2>, unaligned, false >, false >
   , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   typedef SVecSVecCrossExpr<VT1,VT2>  CPE;  //!< Type of the cross product expression.
   typedef typename CPE::ResultType    RT;   //!< Result type of the cross product expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DenseSubvector<CPE,unaligned,false>  This;           //!< Type of this DenseSubvector instance.
   typedef typename SubvectorTrait<RT>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType   TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename CPE::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename CPE::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const ResultType                     CompositeType;  //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DenseSubvector specialization class.
   //
   // \param vector The sparse vector/sparse vector cross product expression.
   // \param index The first index of the subvector in the given expression.
   // \param n The size of the subvector.
   */
   explicit inline DenseSubvector( const CPE& vector, size_t index, size_t n )
      : vector_( vector )  // The sparse vector/sparse vector cross product expression
      , offset_( index  )  // The offset of the subvector within the cross product expression
      , size_  ( n      )  // The size of the subvector
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset_+index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return size_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const {
      return vector_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CPE          vector_;  //!< The sparse vector/sparse vector cross product expression.
   const size_t offset_;  //!< The offset of the subvector within the cross product expression.
   const size_t size_;    //!< The size of the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< bool AF1, typename VT, bool AF2, bool TF >
   friend const DenseSubvector<VT,AF1,TF>
      subvector( const DenseSubvector<VT,AF2,TF>& dv, size_t index, size_t size );

   template< typename VT3, bool AF, bool TF >
   friend bool isIntact( const DenseSubvector<VT3,AF,TF>& dv );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseVector<VT3,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseVector<VT3,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF >
   friend bool isSame( const DenseSubvector<VT3,AF,TF>& a, const DenseSubvector<VT3,AF,TF>& b );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryAddAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool trySubAssign( const DenseSubvector<VT2,AF,TF>& lhs, const Vector<VT3,TF>& rhs, size_t index );

   template< typename VT3, bool AF, bool TF, typename VT4 >
   friend bool tryMultAssign( const DenseSubvector<VT3,AF,TF>& lhs, const Vector<VT4,TF>& rhs, size_t index );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  DENSESUBVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseSubvector operators */
//@{
template< typename VT, bool AF, bool TF >
inline void reset( DenseSubvector<VT,AF,TF>& dv );

template< typename VT, bool AF, bool TF >
inline void clear( DenseSubvector<VT,AF,TF>& dv );

template< typename VT, bool AF, bool TF >
inline bool isDefault( const DenseSubvector<VT,AF,TF>& dv );

template< typename VT, bool AF, bool TF >
inline bool isIntact( const DenseSubvector<VT,AF,TF>& dv );

template< typename VT, bool AF, bool TF >
inline bool isSame( const DenseSubvector<VT,AF,TF>& a, const DenseVector<VT,TF>& b );

template< typename VT, bool AF, bool TF >
inline bool isSame( const DenseVector<VT,TF>& a, const DenseSubvector<VT,AF,TF>& b );

template< typename VT, bool AF, bool TF >
inline bool isSame( const DenseSubvector<VT,AF,TF>& a, const DenseSubvector<VT,AF,TF>& b );
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void reset( DenseSubvector<VT,AF,TF>& dv )
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
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void clear( DenseSubvector<VT,AF,TF>& dv )
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
// the function returns \a true in case all subvector elements are 0 and \a false in case any
// subvector element is not 0. The following example demonstrates the use of the \a isDefault
// function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isDefault( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isDefault( const DenseSubvector<VT,AF,TF>& dv )
{
   for( size_t i=0UL; i<dv.size(); ++i )
      if( !isDefault( dv[i] ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given dense subvector vector are intact.
// \ingroup dense_subvector
//
// \param dv The dense subvector to be tested.
// \return \a true in case the given subvector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dense subvector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isIntact( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isIntact( const DenseSubvector<VT,AF,TF>& dv )
{
   return ( dv.offset_ + dv.size_ <= dv.vector_.size() &&
            isIntact( dv.vector_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense vector and subvector represent the same observable state.
// \ingroup dense_subvector
//
// \param a The dense subvector to be tested for its state.
// \param b The dense vector to be tested for its state.
// \return \a true in case the dense subvector and vector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given subvector refers to the entire range
// of the given dense vector and by that represents the same observable state. In this case, the
// function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isSame( const DenseSubvector<VT,AF,TF>& a, const DenseVector<VT,TF>& b )
{
   return ( isSame( a.vector_, ~b ) && ( a.size() == (~b).size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense vector and subvector represent the same observable state.
// \ingroup dense_subvector
//
// \param a The dense vector to be tested for its state.
// \param b The dense subvector to be tested for its state.
// \return \a true in case the dense vector and subvector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given subvector refers to the entire range
// of the given dense vector and by that represents the same observable state. In this case, the
// function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isSame( const DenseVector<VT,TF>& a, const DenseSubvector<VT,AF,TF>& b )
{
   return ( isSame( ~a, b.vector_ ) && ( (~a).size() == b.size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given subvectors represent the same observable state.
// \ingroup dense_subvector
//
// \param a The first dense subvector to be tested for its state.
// \param b The second dense subvector to be tested for its state.
// \return \a true in case the two subvectors share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given subvectors refer to exactly the
// same range of the same dense vector. In case both subvectors represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isSame( const DenseSubvector<VT,AF,TF>& a, const DenseSubvector<VT,AF,TF>& b )
{
   return ( isSame( a.vector_, b.vector_ ) && ( a.offset_ == b.offset_ ) && ( a.size_ == b.size_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a dense subvector.
// \ingroup dense_subvector
//
// \param lhs The target left-hand side dense subvector.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool tryAssign( const DenseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a dense subvector.
// \ingroup dense_subvector
//
// \param lhs The target left-hand side dense subvector.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool tryAddAssign( const DenseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a dense
//        subvector.
// \ingroup dense_subvector
//
// \param lhs The target left-hand side dense subvector.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool trySubAssign( const DenseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a dense
//        subvector.
// \ingroup dense_subvector
//
// \param lhs The target left-hand side dense subvector.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool tryMultAssign( const DenseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given dense subvector.
// \ingroup dense_subvector
//
// \param dv The subvector to be derestricted.
// \return Subvector without access restrictions.
//
// This function removes all restrictions on the data access to the given subvector. It returns a
// subvector that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DerestrictTrait< DenseSubvector<VT,AF,TF> >::Type
   derestrict( DenseSubvector<VT,AF,TF>& dv )
{
   typedef typename DerestrictTrait< DenseSubvector<VT,AF,TF> >::Type  ReturnType;
   return ReturnType( derestrict( dv.vector_ ), dv.offset_, dv.size_ );
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
template< bool AF1     // Required alignment flag
        , typename VT  // Type of the dense vector
        , bool AF2     // Present alignment flag
        , bool TF >    // Transpose flag
inline const DenseSubvector<VT,AF1,TF>
   subvector( const DenseSubvector<VT,AF2,TF>& dv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   if( index + size > dv.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }

   return DenseSubvector<VT,AF1,TF>( dv.vector_, dv.offset_ + index, size );
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
template< typename VT, bool AF, bool TF >
struct IsRestricted< DenseSubvector<VT,AF,TF> > : public IsTrue< IsRestricted<VT>::value >
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
template< typename VT, bool AF, bool TF >
struct DerestrictTrait< DenseSubvector<VT,AF,TF> >
{
   typedef DenseSubvector< typename RemoveReference< typename DerestrictTrait<VT>::Type >::Type, AF, TF >  Type;
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
template< typename VT, bool AF, bool TF >
struct HasConstDataAccess< DenseSubvector<VT,AF,TF> >
   : public IsTrue< HasConstDataAccess<VT>::value >
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
template< typename VT, bool AF, bool TF >
struct HasMutableDataAccess< DenseSubvector<VT,AF,TF> >
   : public IsTrue< HasMutableDataAccess<VT>::value >
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
template< typename VT, bool TF >
struct IsAligned< DenseSubvector<VT,aligned,TF> > : public IsTrue<true>
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
template< typename VT, bool AF, bool TF, typename T >
struct AddTrait< DenseSubvector<VT,AF,TF>, T >
{
   typedef typename AddTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct AddTrait< T, DenseSubvector<VT,AF,TF> >
{
   typedef typename AddTrait< T, typename SubvectorTrait<VT>::Type >::Type  Type;
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
template< typename VT, bool AF, bool TF, typename T >
struct SubTrait< DenseSubvector<VT,AF,TF>, T >
{
   typedef typename SubTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct SubTrait< T, DenseSubvector<VT,AF,TF> >
{
   typedef typename SubTrait< T, typename SubvectorTrait<VT>::Type >::Type  Type;
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
template< typename VT, bool AF, bool TF, typename T >
struct MultTrait< DenseSubvector<VT,AF,TF>, T >
{
   typedef typename MultTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct MultTrait< T, DenseSubvector<VT,AF,TF> >
{
   typedef typename MultTrait< T, typename SubvectorTrait<VT>::Type >::Type  Type;
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
template< typename VT, bool AF, bool TF, typename T >
struct CrossTrait< DenseSubvector<VT,AF,TF>, T >
{
   typedef typename CrossTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct CrossTrait< T, DenseSubvector<VT,AF,TF> >
{
   typedef typename CrossTrait< T, typename SubvectorTrait<VT>::Type >::Type  Type;
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
template< typename VT, bool AF, bool TF, typename T >
struct DivTrait< DenseSubvector<VT,AF,TF>, T >
{
   typedef typename DivTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct DivTrait< T, DenseSubvector<VT,AF,TF> >
{
   typedef typename DivTrait< T, typename SubvectorTrait<VT>::Type >::Type  Type;
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
template< typename VT, bool AF, bool TF >
struct SubvectorTrait< DenseSubvector<VT,AF,TF> >
{
   typedef typename SubvectorTrait< typename DenseSubvector<VT,AF,TF>::ResultType >::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTOREXPRTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< DenseSubvector<VT,AF1,TF>, AF2 >
{
   typedef DenseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< const DenseSubvector<VT,AF1,TF>, AF2 >
{
   typedef DenseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< volatile DenseSubvector<VT,AF1,TF>, AF2 >
{
   typedef DenseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< const volatile DenseSubvector<VT,AF1,TF>, AF2 >
{
   typedef DenseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool AF >
struct SubvectorExprTrait< DVecDVecCrossExpr<VT1,VT2>, AF >
{
 public:
   //**********************************************************************************************
   typedef DenseSubvector< DVecDVecCrossExpr<VT1,VT2>, unaligned, false >  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool AF >
struct SubvectorExprTrait< DVecSVecCrossExpr<VT1,VT2>, AF >
{
 public:
   //**********************************************************************************************
   typedef DenseSubvector< DVecSVecCrossExpr<VT1,VT2>, unaligned, false >  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool AF >
struct SubvectorExprTrait< SVecDVecCrossExpr<VT1,VT2>, AF >
{
 public:
   //**********************************************************************************************
   typedef DenseSubvector< SVecDVecCrossExpr<VT1,VT2>, unaligned, false >  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool AF >
struct SubvectorExprTrait< SVecSVecCrossExpr<VT1,VT2>, AF >
{
 public:
   //**********************************************************************************************
   typedef DenseSubvector< SVecSVecCrossExpr<VT1,VT2>, unaligned, false >  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
