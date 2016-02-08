//=================================================================================================
/*!
//  \file blaze/math/views/SparseSubvector.h
//  \brief Header file for the SparseSubvector class template
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

#ifndef _BLAZE_MATH_VIEWS_SPARSESUBVECTOR_H_
#define _BLAZE_MATH_VIEWS_SPARSESUBVECTOR_H_


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
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/Subvector.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/util/Assert.h>
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
/*!\defgroup sparse_subvector Sparse Subvector
// \ingroup views
*/
/*!\brief View on a specific subvector of a sparse vector.
// \ingroup sparse_subvector
//
// The SparseSubvector template represents a view on a specific subvector of a sparse vector
// primitive. The type of the sparse vector is specified via the first template parameter:

   \code
   template< typename VT, bool AF, bool TF >
   class SparseSubvector;
   \endcode

//  - VT: specifies the type of the sparse vector primitive. SparseSubvector can be used with
//        every sparse vector primitive or view, but does not work with any vector expression type.
//  - AF: the alignment flag specifies whether the subvector is aligned (\a blaze::aligned) or
//        unaligned (\a blaze::unaligned). The default value is \a blaze::unaligned.
//  - TF: specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//        vector (\a blaze::columnVector). This template parameter doesn't have to be explicitly
//        defined, but is automatically derived from the first template parameter.
//
//
// \n \section sparse_subvector_setup Setup of Sparse Subvectors
//
// A view on a sparse subvector can be created very conveniently via the \c subvector() function:

      \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;

   SparseVectorType x;
   // ... Resizing and initialization

   // Create a subvector from index 8 with a size of 16 (i.e. in the range [8..23])
   blaze::SparseSubvector<SparseVectorType> sv = subvector( x, 8UL, 16UL );
   \endcode

// This view can be treated as any other sparse vector, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The view can also be used on both
// sides of an assignment: The subvector can either be used as an alias to grant write access to
// a specific subvector of a sparse vector primitive on the left-hand side of an assignment or
// to grant read-access to a specific subvector of a sparse vector primitive or expression on
// the right-hand side of an assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>   SparseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   SparseMatrixType A;
   // ... Resizing and initialization

   // Create a subvector from index 0UL with a size of 10 (i.e. in the range [0..9])
   blaze::SparseSubvector<SparseVectorType> sv = subvector( y, 0UL, 10UL );

   // Setting the first ten elements of y to the 2nd row of matrix A
   sv = row( A, 2UL );

   // Setting the second ten elements of y to x
   subvector( y, 10UL, 10UL ) = x;

   // Setting the 3rd row of A to a subvector of y
   row( A, 3UL ) = subvector( y, 3UL, 10UL );

   // Setting y to a subvector of the result of the addition between x and the 1st row of A
   y = subvector( x + row( A, 1UL ), 2UL, 5UL )
   \endcode

// \n \section sparse_subvector_element_access Element access
//
// A sparse subvector can be used like any other sparse vector. For instance, the elements of the
// sparse subvector can be directly accessed with the subscript operator.

   \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  VectorType;
   VectorType v;
   // ... Resizing and initialization

   // Creating an 8-dimensional subvector, starting from index 4
   blaze::SparseSubvector<VectorType> sv = subvector( v, 4UL, 8UL );

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
   typedef blaze::CompressedVector<int,blaze::rowVector>  VectorType;
   typedef blaze::SparseSubvector<VectorType>             SubvectorType;

   VectorType v( 256UL );
   // ... Resizing and initialization

   // Creating a reference to a specific subvector of vector v
   SubvectorType sv = subvector( v, 16UL, 64UL );

   for( SubvectorType::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( SubvectorType::ConstIterator it=sv.begin(); it!=sv.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_subvector_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse subvector can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  VectorType;
   VectorType v( 256UL );  // Non-initialized vector of size 256

   typedef blaze::SparseSubvector<VectorType>  SubvectorType;
   SubvectorType sv( subvector( v, 10UL, 60UL ) );  // View on the range [10..69] of v

   // The subscript operator provides access to all possible elements of the sparse subvector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse subvector, the element is inserted into the
   // subvector.
   sv[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the subvector it is inserted into the subvector, if it is already contained
   // in the subvector its value is modified.
   sv.set( 45UL, -1.2 );

   // An alternative for inserting elements into the subvector is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the subvector.
   sv.insert( 50UL, 3.7 );

   // Just as in case of vectors, elements can also be inserted via the append() function. In
   // case of subvectors, append() also requires that the appended element's index is strictly
   // larger than the currently largest non-zero index of the subvector and that the subvector's
   // capacity is large enough to hold the new element. Note however that due to the nature of
   // a subvector, which may be an alias to the middle of a sparse vector, the append() function
   // does not work as efficiently for a subvector as it does for a vector.
   sv.reserve( 10UL );
   sv.append( 51UL, -2.1 );
   \endcode

// \n \section sparse_subvector_common_operations Common Operations
//
// The current number of subvector elements can be obtained via the \c size() function, the
// current capacity via the \c capacity() function, and the number of non-zero elements via
// the \c nonZeros() function. However, since subvector are views on a specific subvector of
// a vector, several operations are not possible on views, such as resizing and swapping:

   \code
   typedef blaze::CompressedVector<int,blaze::rowVector>  VectorType;
   typedef blaze::SparseSubvector<VectorType>             SubvectorType;

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
// The following example gives an impression of the use of SparseSubvector within arithmetic
// operations. All operations (addition, subtraction, multiplication, scaling, ...) can be
// performed on all possible combinations of dense and sparse vectors with fitting element
// types:

   \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   SparseVectorType s1, s2, s3;
   DenseVectorType d1, d2;

   // ... Resizing and initialization

   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;
   SparseMatrixType A;

   typedef blaze::SparseSubvector<SparseVectorType>  SubvectorType;
   SubvectorType sv( subvector( s1, 0UL, 10UL ) );  // View on the range [0..9] of vector s1

   sv = s2;                           // Sparse vector initialization of the range [0..9]
   subvector( s1, 10UL, 10UL ) = d1;  // Dense vector initialization of the range [10..19]

   s3 = sv + s2;                           // Sparse vector/sparse vector addition
   d2 = d1 + subvector( s1, 10UL, 10UL );  // Dense vector/sparse vector addition
   s2 = sv * subvector( s1, 20UL, 10UL );  // Component-wise vector multiplication

   subvector( s1, 3UL, 4UL ) *= 2.0;      // In-place scaling of the range [3..6]
   s2 = subvector( s1, 7UL, 3UL ) * 2.0;  // Scaling of the range [7..9]
   s2 = 2.0 * subvector( s1, 7UL, 3UL );  // Scaling of the range [7..9]

   subvector( s1, 0UL , 10UL ) += s2;  // Addition assignment
   subvector( s1, 10UL, 10UL ) -= d2;  // Subtraction assignment
   subvector( s1, 20UL, 10UL ) *= sv;  // Multiplication assignment

   double scalar = subvector( s1, 5UL, 10UL ) * trans( d1 );  // Scalar/dot/inner product between two vectors

   A = trans( d1 ) * subvector( s1, 4UL, 16UL );  // Outer product between two vectors
   \endcode

// \n \section sparse_subvector_aligned_subvector Aligned Subvectors
//
// Usually subvectors can be defined anywhere within a vector. They may start at any position and
// may have an arbitrary size (only restricted by the size of the underlying vector). However, in
// contrast to vectors themselves, which are always properly aligned in memory and therefore can
// provide maximum performance, this means that subvectors in general have to be considered to be
// unaligned. This can be made explicit by the \a blaze::unaligned flag:

   \code
   using blaze::unaligned;

   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;

   SparseVectorType x;
   // ... Resizing and initialization

   // Identical creations of an unaligned subvector in the range [8..23]
   blaze::SparseSubvector<SparseVectorType>           sv1 = subvector           ( x, 8UL, 16UL );
   blaze::SparseSubvector<SparseVectorType>           sv2 = subvector<unaligned>( x, 8UL, 16UL );
   blaze::SparseSubvector<SparseVectorType,unaligned> sv3 = subvector           ( x, 8UL, 16UL );
   blaze::SparseSubvector<SparseVectorType,unaligned> sv4 = subvector<unaligned>( x, 8UL, 16UL );
   \endcode

// All of these calls to the \c subvector() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned subvector. Whereas this may provide
// full flexibility in the creation of subvectors, this might result in performance restrictions
// (even in case the specified subvector could be aligned). However, it is also possible to create
// aligned subvectors. Aligned subvectors are identical to unaligned subvectors in all aspects,
// except that they may pose additional alignment restrictions and therefore have less flexibility
// during creation. These restrictions may limit their application, but due to that they don't
// suffer from performance penalties and provide the same performance as the underlying vector.
// Aligned subvectors are created by explicitly specifying the \a blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned subvector in the range [8..23]
   blaze::SparseSubvector<SparseVectorType,aligned> sv = subvector<aligned>( x, 8UL, 16UL );
   \endcode

// In contrast to dense subvectors, which pose several additional alignment restrictions based on
// the used element type, sparse subvectors at this time don't pose any additional restrictions.
// Therefore aligned and unaligned sparse subvectors are truly fully identical. Note however that
// this is not true for dense subvectors (see the DenseSubvector class description)!
//
// \n \section sparse_subvector_on_sparse_subvector Subvectors on Subvectors
//
// It is also possible to create a subvector view on another subvector. In this context it is
// important to remember that the type returned by the \c subvector() function is the same type
// as the type of the given subvector, since the view on a subvector is just another view on the
// underlying sparse vector:

   \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  VectorType;
   typedef blaze::SparseSubvector<VectorType>                SubvectorType;

   VectorType s1;

   // ... Resizing and initialization

   // Creating a subvector view on the sparse vector s1
   SubvectorType sv1 = subvector( s1, 5UL, 10UL );

   // Creating a subvector view on the sparse subvector sv1
   SubvectorType sv2 = subvector( sv1, 1UL, 5UL );
   \endcode
*/
template< typename VT                         // Type of the sparse vector
        , bool AF = unaligned                 // Alignment flag
        , bool TF = IsRowVector<VT>::value >  // Transpose flag
class SparseSubvector : public SparseVector< SparseSubvector<VT,AF,TF>, TF >
                      , private Subvector
{
 private:
   //**Type definitions****************************************************************************
   //! Composite data type of the sparse vector expression.
   typedef typename If< IsExpression<VT>, VT, VT& >::Type  Operand;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SparseSubvector<VT,AF,TF>           This;           //!< Type of this SparseSubvector instance.
   typedef typename SubvectorTrait<VT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename VT::ElementType            ElementType;    //!< Type of the subvector elements.
   typedef typename VT::ReturnType             ReturnType;     //!< Return type for expression template evaluations
   typedef const SparseSubvector&              CompositeType;  //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   typedef typename VT::ConstReference  ConstReference;

   //! Reference to a non-constant subvector value.
   typedef typename If< IsConst<VT>, ConstReference, typename VT::Reference >::Type  Reference;
   //**********************************************************************************************

   //**SubvectorElement class definition***********************************************************
   /*!\brief Access proxy for a specific element of the sparse subvector.
   */
   template< typename VectorType      // Type of the sparse vector
           , typename IteratorType >  // Type of the sparse vector iterator
   class SubvectorElement : private SparseElement
   {
    private:
      //*******************************************************************************************
      //! Compilation switch for the return type of the value member function.
      /*! The \a returnConst compile time constant expression represents a compilation switch for
          the return type of the value member function. In case the given vector type \a VectorType
          is const qualified, \a returnConst will be set to 1 and the value member function will
          return a reference to const. Otherwise \a returnConst will be set to 0 and the value
          member function will offer write access to the sparse vector elements. */
      enum { returnConst = IsConst<VectorType>::value };
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
      /*!\brief Constructor for the SubvectorElement class.
      //
      // \param pos Iterator to the current position within the sparse subvector.
      // \param offset The offset within the according sparse vector.
      */
      inline SubvectorElement( IteratorType pos, size_t offset )
         : pos_   ( pos    )  // Iterator to the current position within the sparse subvector
         , offset_( offset )  // Offset within the according sparse vector
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse subvector element.
      //
      // \param v The new value of the sparse subvector element.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the addition.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the subtraction.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the multiplication.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the division.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse subvector element at the current iterator position.
      //
      // \return Reference to the sparse subvector element at the current iterator position.
      */
      inline const SubvectorElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse subvector element.
      //
      // \return The current value of the sparse subvector element.
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
      IteratorType pos_;  //!< Iterator to the current position within the sparse subvector.
      size_t offset_;     //!< Offset within the according sparse vector.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**SubvectorIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse subvector.
   */
   template< typename VectorType      // Type of the sparse vector
           , typename IteratorType >  // Type of the sparse vector iterator
   class SubvectorIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::forward_iterator_tag                  IteratorCategory;  //!< The iterator category.
      typedef SubvectorElement<VectorType,IteratorType>  ValueType;         //!< Type of the underlying elements.
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
      /*!\brief Default constructor for the SubvectorIterator class.
      */
      inline SubvectorIterator()
         : pos_   ()  // Iterator to the current sparse element
         , offset_()  // The offset of the subvector within the sparse vector
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the SubvectorIterator class.
      //
      // \param iterator Iterator to the current sparse element.
      // \param index The starting index of the subvector within the sparse vector.
      */
      inline SubvectorIterator( IteratorType iterator, size_t index )
         : pos_   ( iterator )  // Iterator to the current sparse element
         , offset_( index    )  // The offset of the subvector within the sparse vector
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubvectorIterator instances.
      //
      // \param it The subvector iterator to be copied.
      */
      template< typename VectorType2, typename IteratorType2 >
      inline SubvectorIterator( const SubvectorIterator<VectorType2,IteratorType2>& it )
         : pos_   ( it.base()   )  // Iterator to the current sparse element.
         , offset_( it.offset() )  // The offset of the subvector within the sparse vector
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubvectorIterator& operator++() {
         ++pos_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubvectorIterator operator++( int ) {
         const SubvectorIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse subvector element.
      //
      // \return Reference to the sparse subvector element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse subvector element.
      //
      // \return Pointer to the sparse subvector element.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side subvector iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename VectorType2, typename IteratorType2 >
      inline bool operator==( const SubvectorIterator<VectorType2,IteratorType2>& rhs ) const {
         return base() == rhs.base();
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side subvector iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename VectorType2, typename IteratorType2 >
      inline bool operator!=( const SubvectorIterator<VectorType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two subvector iterators.
      //
      // \param rhs The right-hand side subvector iterator.
      // \return The number of elements between the two subvector iterators.
      */
      inline DifferenceType operator-( const SubvectorIterator& rhs ) const {
         return pos_ - rhs.pos_;
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the subvector iterator.
      //
      // \return The current position of the subvector iterator.
      */
      inline IteratorType base() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Offset function**************************************************************************
      /*!\brief Access to the offset of the subvector iterator.
      //
      // \return The offset of the subvector iterator.
      */
      inline size_t offset() const {
         return offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;     //!< Iterator to the current sparse element.
      size_t       offset_;  //!< The offset of the subvector within the sparse vector.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   typedef SubvectorIterator<const VT,typename VT::ConstIterator>  ConstIterator;

   //! Iterator over non-constant elements.
   typedef typename If< IsConst<VT>, ConstIterator, SubvectorIterator<VT,typename VT::Iterator> >::Type  Iterator;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = VT::smpAssignable };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SparseSubvector( Operand vector, size_t index, size_t n );
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
                            inline SparseSubvector& operator= ( const SparseSubvector& rhs );
   template< typename VT2 > inline SparseSubvector& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline SparseSubvector& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline SparseSubvector& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline SparseSubvector& operator*=( const Vector<VT2,TF>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseSubvector >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, SparseSubvector >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t           size() const;
                              inline size_t           capacity() const;
                              inline size_t           nonZeros() const;
                              inline void             reset();
                              inline Iterator         set    ( size_t index, const ElementType& value );
                              inline Iterator         insert ( size_t index, const ElementType& value );
                              inline void             erase  ( size_t index );
                              inline Iterator         erase  ( Iterator pos );
                              inline Iterator         erase  ( Iterator first, Iterator last );
                              inline void             reserve( size_t n );
   template< typename Other > inline SparseSubvector& scale  ( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
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

   inline bool canSMPAssign() const;

   template< typename VT2 >   inline void assign   ( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 >   inline void assign   ( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 >   inline void addAssign( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 >   inline void addAssign( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 >   inline void subAssign( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 >   inline void subAssign( const SparseVector<VT2,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      vector_;  //!< The sparse vector containing the subvector.
   const size_t offset_;  //!< The offset of the subvector within the sparse vector.
   const size_t size_;    //!< The size of the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< bool AF1, typename VT2, bool AF2, bool TF2 >
   friend const SparseSubvector<VT2,AF1,TF2>
      subvector( const SparseSubvector<VT2,AF2,TF2>& sv, size_t index, size_t size );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isIntact( const SparseSubvector<VT2,AF2,TF2>& sv );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const SparseSubvector<VT2,AF2,TF2>& a, const SparseVector<VT2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const SparseVector<VT2,TF2>& a, const SparseSubvector<VT2,AF2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2 >
   friend bool isSame( const SparseSubvector<VT2,AF2,TF2>& a, const SparseSubvector<VT2,AF2,TF2>& b );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryAssign( const SparseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryAddAssign( const SparseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool trySubAssign( const SparseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2, typename VT3 >
   friend bool tryMultAssign( const SparseSubvector<VT2,AF2,TF2>& lhs, const Vector<VT3,TF2>& rhs, size_t index );

   template< typename VT2, bool AF2, bool TF2 >
   friend typename DerestrictTrait< SparseSubvector<VT2,AF2,TF2> >::Type
      derestrict( SparseSubvector<VT2,AF2,TF2>& sv );
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE  ( VT );
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
/*!\brief The constructor for SparseSubvector.
//
// \param vector The sparse vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param n The size of the subvector.
// \exception std::invalid_argument Invalid subvector specification.
//
// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the size of the given vector or the subvector is specified beyond the size of the vector)
// a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline SparseSubvector<VT,AF,TF>::SparseSubvector( Operand vector, size_t index, size_t n )
   : vector_( vector )  // The sparse vector containing the subvector
   , offset_( index  )  // The offset of the subvector within the sparse vector
   , size_  ( n      )  // The size of the subvector
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
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Reference
   SparseSubvector<VT,AF,TF>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return vector_[offset_+index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstReference
   SparseSubvector<VT,AF,TF>::operator[]( size_t index ) const
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
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Reference
   SparseSubvector<VT,AF,TF>::at( size_t index )
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
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstReference
   SparseSubvector<VT,AF,TF>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator SparseSubvector<VT,AF,TF>::begin()
{
   if( offset_ == 0UL )
      return Iterator( vector_.begin(), offset_ );
   else
      return Iterator( vector_.lowerBound( offset_ ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator SparseSubvector<VT,AF,TF>::begin() const
{
   if( offset_ == 0UL )
      return ConstIterator( vector_.cbegin(), offset_ );
   else
      return ConstIterator( vector_.lowerBound( offset_ ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator SparseSubvector<VT,AF,TF>::cbegin() const
{
   if( offset_ == 0UL )
      return ConstIterator( vector_.cbegin(), offset_ );
   else
      return ConstIterator( vector_.lowerBound( offset_ ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator SparseSubvector<VT,AF,TF>::end()
{
   if( offset_ + size_ == vector_.size() )
      return Iterator( vector_.end(), offset_ );
   else
      return Iterator( vector_.lowerBound( offset_ + size_ ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator SparseSubvector<VT,AF,TF>::end() const
{
   if( offset_ + size_ == vector_.size() )
      return ConstIterator( vector_.cend(), offset_ );
   else
      return ConstIterator( vector_.lowerBound( offset_ + size_ ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator SparseSubvector<VT,AF,TF>::cend() const
{
   if( offset_ + size_ == vector_.size() )
      return ConstIterator( vector_.cend(), offset_ );
   else
      return ConstIterator( vector_.lowerBound( offset_ + size_ ), offset_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for SparseSubvector.
//
// \param rhs Sparse subvector to be copied.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Subvector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline SparseSubvector<VT,AF,TF>&
   SparseSubvector<VT,AF,TF>::operator=( const SparseSubvector& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &vector_ == &rhs.vector_ && offset_ == rhs.offset_ ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !tryAssign( vector_, rhs, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   if( rhs.canAlias( &vector_ ) ) {
      const ResultType tmp( rhs );
      reset();
      assign( left, tmp );
   }
   else {
      reset();
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Dense vector to be assigned.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline SparseSubvector<VT,AF,TF>&
   SparseSubvector<VT,AF,TF>::operator=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

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

   if( IsReference<Right>::value || right.canAlias( &vector_ ) ) {
      const typename VT2::ResultType tmp( right );
      reset();
      assign( left, tmp );
   }
   else {
      reset();
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the sparse subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline SparseSubvector<VT,AF,TF>&
   SparseSubvector<VT,AF,TF>::operator+=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   typedef typename AddTrait<ResultType,typename VT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (~rhs) );

   if( !tryAssign( vector_, tmp, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the sparse subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline SparseSubvector<VT,AF,TF>&
   SparseSubvector<VT,AF,TF>::operator-=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   typedef typename SubTrait<ResultType,typename VT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (~rhs) );

   if( !tryAssign( vector_, tmp, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side vector
inline SparseSubvector<VT,AF,TF>&
   SparseSubvector<VT,AF,TF>::operator*=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( typename VT2::ResultType );

   typedef typename MultTrait<ResultType,typename VT2::ResultType>::Type  MultType;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( MultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (~rhs) );

   if( !tryAssign( vector_, tmp, offset_ ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   typename DerestrictTrait<This>::Type left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a sparse subvector
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the assigned subvector.
//
// This operator can only be used for built-in data types. Additionally, the elements of
// the sparse subvector must support the multiplication assignment operator for the given
// scalar built-in data type.
*/
template< typename VT       // Type of the sparse vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubvector<VT,AF,TF> >::Type&
   SparseSubvector<VT,AF,TF>::operator*=( Other rhs )
{
   const Iterator last( end() );
   for( Iterator element=begin(); element!=last; ++element )
      element->value() *= rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse subvector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the assigned subvector.
//
// This operator can only be used for built-in data types. Additionally, the elements of the
// sparse subvector must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
*/
template< typename VT       // Type of the sparse vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, SparseSubvector<VT,AF,TF> >::Type&
   SparseSubvector<VT,AF,TF>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename DivTrait<ElementType,Other>::Type  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   const Iterator last( end() );

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( Iterator element=begin(); element!=last; ++element )
         element->value() *= tmp;
   }
   else {
      for( Iterator element=begin(); element!=last; ++element )
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
/*!\brief Returns the size/dimension of the sparse subvector.
//
// \return The size of the sparse subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline size_t SparseSubvector<VT,AF,TF>::size() const
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the sparse subvector.
//
// \return The capacity of the sparse subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline size_t SparseSubvector<VT,AF,TF>::capacity() const
{
   return nonZeros() + vector_.capacity() - vector_.nonZeros();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the subvector.
//
// \return The number of non-zero elements in the subvector.
//
// Note that the number of non-zero elements is always smaller than the size of the subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline size_t SparseSubvector<VT,AF,TF>::nonZeros() const
{
   return end() - begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void SparseSubvector<VT,AF,TF>::reset()
{
   vector_.erase( vector_.lowerBound( offset_ ), vector_.lowerBound( offset_ + size_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an element of the sparse subvector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse subvector. In case the sparse subvector
// already contains an element with index \a index its value is modified, else a new element with
// the given \a value is inserted.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator
   SparseSubvector<VT,AF,TF>::set( size_t index, const ElementType& value )
{
   return Iterator( vector_.set( offset_ + index, value ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the sparse subvector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse subvector access index.
//
// This function inserts a new element into the sparse subvector. However, duplicate elements
// are not allowed. In case the sparse subvector already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator
   SparseSubvector<VT,AF,TF>::insert( size_t index, const ElementType& value )
{
   return Iterator( vector_.insert( offset_ + index, value ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse subvector.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void SparseSubvector<VT,AF,TF>::erase( size_t index )
{
   vector_.erase( offset_ + index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the sparse subvector.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator SparseSubvector<VT,AF,TF>::erase( Iterator pos )
{
   return Iterator( vector_.erase( pos.base() ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the sparse subvector.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse subvector.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator
   SparseSubvector<VT,AF,TF>::erase( Iterator first, Iterator last )
{
   return Iterator( vector_.erase( first.base(), last.base() ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the sparse subvector.
//
// \param n The new minimum capacity of the sparse subvector.
// \return void
//
// This function increases the capacity of the sparse subvector to at least \a n elements. The
// current values of the subvector elements are preserved.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
void SparseSubvector<VT,AF,TF>::reserve( size_t n )
{
   const size_t current( capacity() );

   if( n > current ) {
      vector_.reserve( vector_.capacity() + n - current );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the sparse subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the sparse subvector.
*/
template< typename VT       // Type of the sparse vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the scalar value
inline SparseSubvector<VT,AF,TF>& SparseSubvector<VT,AF,TF>::scale( const Other& scalar )
{
   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific subvector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// subvector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse subvector (the end() iterator) is returned. Note that
// the returned sparse subvector iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator
   SparseSubvector<VT,AF,TF>::find( size_t index )
{
   const typename VT::Iterator pos( vector_.find( offset_ + index ) );

   if( pos != vector_.end() )
      return Iterator( pos, offset_ );
   else
      return end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific subvector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// subvector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse subvector (the end() iterator) is returned. Note that
// the returned sparse subvector iterator is subject to invalidation due to inserting operations
// via the subscript operator or the insert() function!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator
   SparseSubvector<VT,AF,TF>::find( size_t index ) const
{
   const typename VT::ConstIterator pos( vector_.find( offset_ + index ) );

   if( pos != vector_.end() )
      return Iterator( pos, offset_ );
   else
      return end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator
// or the insert() function!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator
   SparseSubvector<VT,AF,TF>::lowerBound( size_t index )
{
   return Iterator( vector_.lowerBound( offset_ + index ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator
// or the insert() function!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator
   SparseSubvector<VT,AF,TF>::lowerBound( size_t index ) const
{
   return ConstIterator( vector_.lowerBound( offset_ + index ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator
// or the insert() function!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::Iterator
   SparseSubvector<VT,AF,TF>::upperBound( size_t index )
{
   return Iterator( vector_.upperBound( offset_ + index ), offset_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator
// or the insert() function!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename SparseSubvector<VT,AF,TF>::ConstIterator
   SparseSubvector<VT,AF,TF>::upperBound( size_t index ) const
{
   return ConstIterator( vector_.upperBound( offset_ + index ), offset_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  LOW-LEVEL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Appending an element to the sparse subvector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse subvector with elements. It
// appends a new element to the end of the sparse subvector without any memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse subvector
//  - the current number of non-zero elements must be smaller than the capacity of the subvector
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void SparseSubvector<VT,AF,TF>::append( size_t index, const ElementType& value, bool check )
{
   if( offset_ + size_ == vector_.size() )
      vector_.append( offset_ + index, value, check );
   else if( !check || !isDefault( value ) )
      vector_.insert( offset_ + index, value );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the sparse subvector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse subvector, \a false if not.
//
// This function returns whether the given address can alias with the sparse subvector. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename VT       // Type of the sparse vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubvector<VT,AF,TF>::canAlias( const Other* alias ) const
{
   return vector_.isAliased( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the sparse subvector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse subvector, \a false if not.
//
// This function returns whether the given address is aliased with the sparse subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the sparse vector
        , bool AF           // Alignment flag
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool SparseSubvector<VT,AF,TF>::isAliased( const Other* alias ) const
{
   return vector_.isAliased( alias );
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
// vector).
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool SparseSubvector<VT,AF,TF>::canSMPAssign() const
{
   return false;
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
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void SparseSubvector<VT,AF,TF>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   reserve( (~rhs).size() );

   for( size_t i=0UL; i<size(); ++i ) {
      append( i, (~rhs)[i], true );
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
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side sparse vector
inline void SparseSubvector<VT,AF,TF>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   reserve( (~rhs).nonZeros() );

   for( typename VT2::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element ) {
      append( element->index(), element->value(), true );
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
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void SparseSubvector<VT,AF,TF>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   reset();
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
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side sparse vector
inline void SparseSubvector<VT,AF,TF>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   typedef typename AddTrait<ResultType,typename VT2::ResultType>::Type  AddType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (~rhs) ) );
   reset();
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
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side dense vector
inline void SparseSubvector<VT,AF,TF>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   reset();
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
template< typename VT     // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename VT2 >  // Type of the right-hand side sparse vector
inline void SparseSubvector<VT,AF,TF>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   typedef typename SubTrait<ResultType,typename VT2::ResultType>::Type  SubType;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (~rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (~rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************








//=================================================================================================
//
//  SPARSESUBVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseSubvector operators */
//@{
template< typename VT, bool AF, bool TF >
inline void reset( SparseSubvector<VT,AF,TF>& sv );

template< typename VT, bool AF, bool TF >
inline void clear( SparseSubvector<VT,AF,TF>& sv );

template< typename VT, bool AF, bool TF >
inline bool isDefault( const SparseSubvector<VT,AF,TF>& sv );

template< typename VT, bool AF, bool TF >
inline bool isIntact( const SparseSubvector<VT,AF,TF>& sv );

template< typename VT, bool AF, bool TF >
inline bool isSame( const SparseSubvector<VT,AF,TF>& a, const SparseVector<VT,TF>& b );

template< typename VT, bool AF, bool TF >
inline bool isSame( const SparseVector<VT,TF>& a, const SparseSubvector<VT,AF,TF>& b );

template< typename VT, bool AF, bool TF >
inline bool isSame( const SparseSubvector<VT,AF,TF>& a, const SparseSubvector<VT,AF,TF>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given sparse subvector.
// \ingroup sparse_subvector
//
// \param sv The sparse subvector to be resetted.
// \return void
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void reset( SparseSubvector<VT,AF,TF>& sv )
{
   sv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given sparse subvector.
// \ingroup sparse_subvector
//
// \param sv The sparse subvector to be cleared.
// \return void
//
// Clearing a sparse subvector is equivalent to resetting it via the reset() function.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void clear( SparseSubvector<VT,AF,TF>& sv )
{
   sv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse subvector is in default state.
// \ingroup sparse_subvector
//
// \param sv The sparse subvector to be tested for its default state.
// \return \a true in case the given subvector is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse subvector is in default state. For instance, in case
// the subvector is instantiated for a vector of built-in integral or floating point data type,
// the function returns \a true in case all subvector elements are 0 and \a false in case any
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedVector<double,rowVector> v;
   // ... Resizing and initialization
   if( isDefault( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isDefault( const SparseSubvector<VT,AF,TF>& sv )
{
   typedef typename SparseSubvector<VT,AF,TF>::ConstIterator  ConstIterator;

   const ConstIterator end( sv.end() );
   for( ConstIterator element=sv.begin(); element!=end; ++element )
      if( !isDefault( element->value() ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given sparse subvector vector are intact.
// \ingroup sparse_subvector
//
// \param sv The sparse subvector to be tested.
// \return \a true in case the given subvector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the sparse subvector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::CompressedVector<double,rowVector> v;
   // ... Resizing and initialization
   if( isIntact( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isIntact( const SparseSubvector<VT,AF,TF>& sv )
{
   return ( sv.offset_ + sv.size_ <= sv.vector_.size() &&
            isIntact( sv.vector_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse vector and subvector represent the same observable state.
// \ingroup sparse_subvector
//
// \param a The sparse subvector to be tested for its state.
// \param b The sparse vector to be tested for its state.
// \return \a true in case the sparse subvector and vector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given subvector refers to the entire range
// of the given sparse vector and by that represents the same observable state. In this case, the
// function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isSame( const SparseSubvector<VT,AF,TF>& a, const SparseVector<VT,TF>& b )
{
   return ( isSame( a.vector_, ~b ) && ( a.size() == (~b).size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given sparse vector and subvector represent the same observable state.
// \ingroup sparse_subvector
//
// \param a The sparse vector to be tested for its state.
// \param b The sparse subvector to be tested for its state.
// \return \a true in case the sparse vector and subvector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given subvector refers to the entire range
// of the given sparse vector and by that represents the same observable state. In this case, the
// function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isSame( const SparseVector<VT,TF>& a, const SparseSubvector<VT,AF,TF>& b )
{
   return ( isSame( ~a, b.vector_ ) && ( (~a).size() == b.size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given subvectors represent the same observable state.
// \ingroup sparse_subvector
//
// \param a The first sparse subvector to be tested for its state.
// \param b The second sparse subvector to be tested for its state.
// \return \a true in case the two subvectors share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given subvectors refer to exactly the
// same range of the same sparse vector. In case both subvectors represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline bool isSame( const SparseSubvector<VT,AF,TF>& a, const SparseSubvector<VT,AF,TF>& b )
{
   return ( isSame( a.vector_, b.vector_ ) && ( a.offset_ == b.offset_ ) && ( a.size_ == b.size_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a sparse subvector.
// \ingroup sparse_subvector
//
// \param lhs The target left-hand side sparse subvector.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool tryAssign( const SparseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a sparse subvector.
// \ingroup sparse_subvector
//
// \param lhs The target left-hand side sparse subvector.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool tryAddAssign( const SparseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a sparse
//        subvector.
// \ingroup sparse_subvector
//
// \param lhs The target left-hand side sparse subvector.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool trySubAssign( const SparseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a sparse
//        subvector.
// \ingroup sparse_subvector
//
// \param lhs The target left-hand side sparse subvector.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the sparse vector
        , bool AF         // Alignment flag
        , bool TF         // Transpose flag
        , typename VT2 >  // Type of the right-hand side vector
inline bool tryMultAssign( const SparseSubvector<VT1,AF,TF>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.vector_, ~rhs, lhs.offset_ + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given sparse subvector.
// \ingroup sparse_subvector
//
// \param sv The subvector to be derestricted.
// \return Subvector without access restrictions.
//
// This function removes all restrictions on the data access to the given subvector. It returns a
// subvector that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the sparse vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline typename DerestrictTrait< SparseSubvector<VT,AF,TF> >::Type
   derestrict( SparseSubvector<VT,AF,TF>& sv )
{
   typedef typename DerestrictTrait< SparseSubvector<VT,AF,TF> >::Type  ReturnType;
   return ReturnType( derestrict( sv.vector_ ), sv.offset_, sv.size_ );
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
/*!\brief Creating a view on a specific subvector of another sparse subvector.
// \ingroup views
//
// \param sv The constant sparse subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the other sparse subvector.
//
// This function returns an expression representing the specified subvector of the given
// sparse subvector.
*/
template< bool AF1     // Required alignment flag
        , typename VT  // Type of the sparse vector
        , bool AF2     // Present alignment flag
        , bool TF >    // Transpose flag
inline const SparseSubvector<VT,AF1,TF>
   subvector( const SparseSubvector<VT,AF2,TF>& sv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   if( index + size > sv.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }

   return SparseSubvector<VT,AF1,TF>( sv.vector_, sv.offset_ + index, size );
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
struct IsRestricted< SparseSubvector<VT,AF,TF> > : public If< IsRestricted<VT>, TrueType, FalseType >::Type
{
   enum { value = IsRestricted<VT>::value };
   typedef typename If< IsRestricted<VT>, TrueType, FalseType >::Type  Type;
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
template< typename VT, bool AF, bool TF >
struct DerestrictTrait< SparseSubvector<VT,AF,TF> >
{
   typedef SparseSubvector< typename RemoveReference< typename DerestrictTrait<VT>::Type >::Type, AF, TF >  Type;
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
template< typename VT, bool AF, bool TF, typename T >
struct AddTrait< SparseSubvector<VT,AF,TF>, T >
{
   typedef typename AddTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct AddTrait< T, SparseSubvector<VT,AF,TF> >
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
struct SubTrait< SparseSubvector<VT,AF,TF>, T >
{
   typedef typename SubTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct SubTrait< T, SparseSubvector<VT,AF,TF> >
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
struct MultTrait< SparseSubvector<VT,AF,TF>, T >
{
   typedef typename MultTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct MultTrait< T, SparseSubvector<VT,AF,TF> >
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
struct CrossTrait< SparseSubvector<VT,AF,TF>, T >
{
   typedef typename CrossTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct CrossTrait< T, SparseSubvector<VT,AF,TF> >
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
struct DivTrait< SparseSubvector<VT,AF,TF>, T >
{
   typedef typename DivTrait< typename SubvectorTrait<VT>::Type, T >::Type  Type;
};

template< typename T, typename VT, bool AF, bool TF >
struct DivTrait< T, SparseSubvector<VT,AF,TF> >
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
struct SubvectorTrait< SparseSubvector<VT,AF,TF> >
{
   typedef typename SubvectorTrait< typename SparseSubvector<VT,AF,TF>::ResultType >::Type  Type;
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
struct SubvectorExprTrait< SparseSubvector<VT,AF1,TF>, AF2 >
{
   typedef SparseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< const SparseSubvector<VT,AF1,TF>, AF2 >
{
   typedef SparseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< volatile SparseSubvector<VT,AF1,TF>, AF2 >
{
   typedef SparseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF1, bool TF, bool AF2 >
struct SubvectorExprTrait< const volatile SparseSubvector<VT,AF1,TF>, AF2 >
{
   typedef SparseSubvector<VT,AF2,TF>  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
