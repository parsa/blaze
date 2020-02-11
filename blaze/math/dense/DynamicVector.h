//=================================================================================================
/*!
//  \file blaze/math/dense/DynamicVector.h
//  \brief Header file for the implementation of an arbitrarily sized vector
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_DENSE_DYNAMICVECTOR_H_
#define _BLAZE_MATH_DENSE_DYNAMICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include <algorithm>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/DenseIterator.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SolveTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Restrict.h>
#include <blaze/system/Thresholds.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/algorithms/Transfer.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/typetraits/RemoveConst.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dynamic_vector DynamicVector
// \ingroup dense_vector
*/
/*!\brief Efficient implementation of an arbitrary sized vector.
// \ingroup dynamic_vector
//
// The DynamicVector class template is the representation of an arbitrary sized vector with
// dynamically allocated elements of arbitrary type. The type of the elements and the transpose
// flag of the vector can be specified via the two template parameters:

   \code
   template< typename Type, bool TF >
   class DynamicVector;
   \endcode

//  - Type: specifies the type of the vector elements. DynamicVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//          vector (\a blaze::columnVector). The default value is \a blaze::columnVector.
//
// These contiguously stored elements can be directly accessed with the subscript operator. The
// numbering of the vector elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right)\f]

// The use of DynamicVector is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// vectors with fitting element types. The following example gives an impression of the use of
// DynamicVector:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::DynamicMatrix;

   DynamicVector<double> a( 2 );  // Non-initialized 2D vector of size 2
   a[0] = 1.0;                    // Initialization of the first element
   a[1] = 2.0;                    // Initialization of the second element

   DynamicVector<double>   b( 2, 2.0  );  // Directly, homogeneously initialized 2D vector
   CompressedVector<float> c( 2 );        // Empty sparse single precision vector
   DynamicVector<double>   d;             // Default constructed dynamic vector
   DynamicMatrix<double>   A;             // Default constructed row-major matrix

   d = a + b;  // Vector addition between vectors of equal element type
   d = a - c;  // Vector subtraction between a dense and sparse vector with different element types
   d = a * b;  // Component-wise vector multiplication

   a *= 2.0;      // In-place scaling of vector
   d  = a * 2.0;  // Scaling of vector a
   d  = 2.0 * a;  // Scaling of vector a

   d += a - b;  // Addition assignment
   d -= a + c;  // Subtraction assignment
   d *= a * b;  // Multiplication assignment

   double scalar = trans( a ) * b;  // Scalar/dot/inner product between two vectors

   A = a * trans( b );  // Outer product between two vectors
   \endcode
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
class DynamicVector
   : public DenseVector< DynamicVector<Type,TF>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This          = DynamicVector<Type,TF>;    //!< Type of this DynamicVector instance.
   using BaseType      = DenseVector<This,TF>;      //!< Base type of this DynamicVector instance.
   using ResultType    = This;                      //!< Result type for expression template evaluations.
   using TransposeType = DynamicVector<Type,!TF>;   //!< Transpose type for expression template evaluations.
   using ElementType   = Type;                      //!< Type of the vector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the vector elements.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations
   using CompositeType = const DynamicVector&;      //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;  //!< Reference to a constant vector value.
   using Pointer        = Type*;        //!< Pointer to a non-constant vector value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant vector value.

   using Iterator      = DenseIterator<Type,aligned>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,aligned>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a DynamicVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using Other = DynamicVector<NewType,TF>;  //!< The type of the other DynamicVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a DynamicVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using Other = DynamicVector<Type,TF>;  //!< The type of the other DynamicVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the vector is involved
       in can be optimized via SIMD operationss. In case the element type of the vector is a
       vectorizable data type, the \a simdEnabled compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool simdEnabled = IsVectorizable_v<Type>;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
            inline DynamicVector() noexcept;
   explicit inline DynamicVector( size_t n );
            inline DynamicVector( size_t n, const Type& init );
            inline DynamicVector( initializer_list<Type> list );

   template< typename Other >
   inline DynamicVector( size_t n, const Other* array );

   template< typename Other, size_t Dim >
   inline DynamicVector( const Other (&array)[Dim] );

   template< typename Other, size_t Dim >
   inline DynamicVector( const std::array<Other,Dim>& array );

                           inline DynamicVector( const DynamicVector& v );
                           inline DynamicVector( DynamicVector&& v ) noexcept;
   template< typename VT > inline DynamicVector( const Vector<VT,TF>& v );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~DynamicVector();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index ) noexcept;
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin () noexcept;
   inline ConstIterator  begin () const noexcept;
   inline ConstIterator  cbegin() const noexcept;
   inline Iterator       end   () noexcept;
   inline ConstIterator  end   () const noexcept;
   inline ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline DynamicVector& operator=( const Type& rhs );
   inline DynamicVector& operator=( initializer_list<Type> list );

   template< typename Other, size_t Dim >
   inline DynamicVector& operator=( const Other (&array)[Dim] );

   template< typename Other, size_t Dim >
   inline DynamicVector& operator=( const std::array<Other,Dim>& array );

   inline DynamicVector& operator=( const DynamicVector& rhs );
   inline DynamicVector& operator=( DynamicVector&& rhs ) noexcept;

   template< typename VT > inline DynamicVector& operator= ( const Vector<VT,TF>& rhs );
   template< typename VT > inline DynamicVector& operator+=( const Vector<VT,TF>& rhs );
   template< typename VT > inline DynamicVector& operator-=( const Vector<VT,TF>& rhs );
   template< typename VT > inline DynamicVector& operator*=( const Vector<VT,TF>& rhs );
   template< typename VT > inline DynamicVector& operator/=( const DenseVector<VT,TF>& rhs );
   template< typename VT > inline DynamicVector& operator%=( const Vector<VT,TF>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   clear();
   inline void   resize( size_t n, bool preserve=true );
   inline void   extend( size_t n, bool preserve=true );
   inline void   reserve( size_t n );
   inline void   shrinkToFit();
   inline void   swap( DynamicVector& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline DynamicVector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT::simdEnabled &&
        IsSIMDCombinable_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDAdd_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDSub_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDMult_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDDiv_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Debugging functions*************************************************************************
   /*!\name Debugging functions */
   //@{
   inline bool isIntact() const noexcept;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT >
   inline auto assign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t addPadding( size_t value ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;             //!< The current size/dimension of the vector.
   size_t capacity_;         //!< The maximum capacity of the vector.
   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated vector elements.
                             /*!< Access to the vector elements is gained via the subscript operator.
                                  The order of the elements is
                                  \f[\left(\begin{array}{*{5}{c}}
                                  0 & 1 & 2 & \cdots & N-1 \\
                                  \end{array}\right)\f] */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEDUCTION GUIDES
//
//=================================================================================================

//*************************************************************************************************
#if BLAZE_CPP17_MODE

template< typename Type >
DynamicVector( size_t, Type* ) -> DynamicVector< RemoveCV_t<Type> >;

template< typename Type, size_t N >
DynamicVector( Type (&)[N] ) -> DynamicVector< RemoveCV_t<Type> >;

template< typename Type, size_t N >
DynamicVector( std::array<Type,N> ) -> DynamicVector<Type>;

#endif
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for DynamicVector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::DynamicVector() noexcept
   : size_    ( 0UL )      // The current size/dimension of the vector
   , capacity_( 0UL )      // The maximum capacity of the vector
   , v_       ( nullptr )  // The vector elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a vector of size \a n. No element initialization is performed!
//
// \param n The size of the vector.
//
// \note This constructor is only responsible to allocate the required dynamic memory. No
// element initialization is performed!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::DynamicVector( size_t n )
   : size_    ( n )                            // The current size/dimension of the vector
   , capacity_( addPadding( n ) )              // The maximum capacity of the vector
   , v_       ( allocate<Type>( capacity_ ) )  // The vector elements
{
   if( IsVectorizable_v<Type> ) {
      for( size_t i=size_; i<capacity_; ++i )
         v_[i] = Type();
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogeneous initialization of all \a n vector elements.
//
// \param n The size of the vector.
// \param init The initial value of the vector elements.
//
// All vector elements are initialized with the specified value.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::DynamicVector( size_t n, const Type& init )
   : DynamicVector( n )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] = init;

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all vector elements.
//
// \param list The initializer list.
//
// This constructor provides the option to explicitly initialize the elements of the vector
// within a constructor call:

   \code
   blaze::DynamicVector<double> v1{ 4.2, 6.3, -1.2 };
   \endcode

// The vector is sized according to the size of the initializer list and all its elements are
// (copy) assigned the elements of the given initializer list.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::DynamicVector( initializer_list<Type> list )
   : DynamicVector( list.size() )
{
   std::fill( std::copy( list.begin(), list.end(), begin() ), end(), Type() );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all vector elements.
//
// \param n The size of the vector.
// \param array Dynamic array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the vector with a
// dynamic array:

   \code
   double* array = new double[4];
   // ... Initialization of the dynamic array
   blaze::DynamicVector<double> v( array, 4UL );
   delete[] array;
   \endcode

// The vector is sized according to the specified size of the array and initialized with the
// values from the given array. Note that it is expected that the given \a array has at least
// \a n elements. Providing an array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the initialization array
inline DynamicVector<Type,TF>::DynamicVector( size_t n, const Other* array )
   : DynamicVector( n )
{
   for( size_t i=0UL; i<n; ++i )
      v_[i] = array[i];

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all vector elements.
//
// \param array Static array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the vector with a
// static array:

   \code
   const int init[4] = { 1, 2, 3 };
   blaze::DynamicVector<int> v( init );
   \endcode

// The vector is sized according to the size of the static array and initialized with the values
// from the given static array. Missing values are initialized with default values (as e.g. the
// fourth element in the example).
*/
template< typename Type   // Data type of the vector
        , bool TF >       // Transpose flag
template< typename Other  // Data type of the static array
        , size_t Dim >    // Dimension of the static array
inline DynamicVector<Type,TF>::DynamicVector( const Other (&array)[Dim] )
   : DynamicVector( Dim )
{
   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of all vector elements from the given std::array.
//
// \param array The given std::array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the vector with a
// std::array:

   \code
   const std::array<int,4UL> init{ 1, 2, 3 };
   blaze::DynamicVector<int> v( init );
   \endcode

// The vector is sized according to the size of the std::array and initialized with the values
// from the given std::array. Missing values are initialized with default values (as e.g. the
// fourth element in the example).
*/
template< typename Type   // Data type of the vector
        , bool TF >       // Transpose flag
template< typename Other  // Data type of the std::array
        , size_t Dim >    // Dimension of the std::array
inline DynamicVector<Type,TF>::DynamicVector( const std::array<Other,Dim>& array )
   : DynamicVector( Dim )
{
   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for DynamicVector.
//
// \param v Vector to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::DynamicVector( const DynamicVector& v )
   : DynamicVector( v.size_ )
{
   BLAZE_INTERNAL_ASSERT( capacity_ <= v.capacity_, "Invalid capacity estimation" );

   smpAssign( *this, ~v );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for DynamicVector.
//
// \param v The vector to be moved into this instance.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::DynamicVector( DynamicVector&& v ) noexcept
   : size_    ( v.size_     )  // The current size/dimension of the vector
   , capacity_( v.capacity_ )  // The maximum capacity of the vector
   , v_       ( v.v_        )  // The vector elements
{
   v.size_     = 0UL;
   v.capacity_ = 0UL;
   v.v_        = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different vectors.
//
// \param v Vector to be copied.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the foreign vector
inline DynamicVector<Type,TF>::DynamicVector( const Vector<VT,TF>& v )
   : DynamicVector( (~v).size() )
{
   if( IsSparseVector_v<VT> ) {
      for( size_t i=0UL; i<size_; ++i ) {
         v_[i] = Type();
      }
   }

   smpAssign( *this, ~v );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for DynamicVector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>::~DynamicVector()
{
   deallocate( v_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::Reference
   DynamicVector<Type,TF>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstReference
   DynamicVector<Type,TF>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::Reference
   DynamicVector<Type,TF>::at( size_t index )
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstReference
   DynamicVector<Type,TF>::at( size_t index ) const
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::Pointer DynamicVector<Type,TF>::data() noexcept
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstPointer DynamicVector<Type,TF>::data() const noexcept
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the dynamic vector.
//
// \return Iterator to the first element of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::Iterator DynamicVector<Type,TF>::begin() noexcept
{
   return Iterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the dynamic vector.
//
// \return Iterator to the first element of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstIterator DynamicVector<Type,TF>::begin() const noexcept
{
   return ConstIterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the dynamic vector.
//
// \return Iterator to the first element of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstIterator DynamicVector<Type,TF>::cbegin() const noexcept
{
   return ConstIterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the dynamic vector.
//
// \return Iterator just past the last element of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::Iterator DynamicVector<Type,TF>::end() noexcept
{
   return Iterator( v_ + size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the dynamic vector.
//
// \return Iterator just past the last element of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstIterator DynamicVector<Type,TF>::end() const noexcept
{
   return ConstIterator( v_ + size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the dynamic vector.
//
// \return Iterator just past the last element of the dynamic vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename DynamicVector<Type,TF>::ConstIterator DynamicVector<Type,TF>::cend() const noexcept
{
   return ConstIterator( v_ + size_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Homogenous assignment to all vector elements.
//
// \param rhs Scalar value to be assigned to all vector elements.
// \return Reference to the assigned vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( const Type& rhs )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] = rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List assignment to all vector elements.
//
// \param list The initializer list.
//
// This assignment operator offers the option to directly assign to all elements of the vector
// by means of an initializer list:

   \code
   blaze::DynamicVector<double> v;
   v = { 4.2, 6.3, -1.2 };
   \endcode

// The vector is resized according to the size of the initializer list and all its elements are
// (copy) assigned the values from the given initializer list.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( initializer_list<Type> list )
{
   resize( list.size(), false );
   std::copy( list.begin(), list.end(), v_ );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all vector elements.
//
// \param array Static array for the assignment.
// \return Reference to the assigned vector.
//
// This assignment operator offers the option to directly set all elements of the vector:

   \code
   const int init[4] = { 1, 2, 3 };
   blaze::DynamicVector<int> v;
   v = init;
   \endcode

// The vector is resized according to the size of the static array and assigned the values from
// the given static array. Missing values are initialized with default values (as e.g. the fourth
// element in the example).
*/
template< typename Type   // Data type of the vector
        , bool TF >       // Transpose flag
template< typename Other  // Data type of the static array
        , size_t Dim >    // Dimension of the static array
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( const Other (&array)[Dim] )
{
   resize( Dim, false );

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all vector elements.
//
// \param array The given std::array for the assignment.
// \return Reference to the assigned vector.
//
// This assignment operator offers the option to directly set all elements of the vector:

   \code
   const std::array<int,4UL> init{ 1, 2, 3 };
   blaze::DynamicVector<int> v;
   v = init;
   \endcode

// The vector is resized according to the size of the std::array and assigned the values from
// the given std::array. Missing values are initialized with default values (as e.g. the fourth
// element in the example).
*/
template< typename Type   // Data type of the vector
        , bool TF >       // Transpose flag
template< typename Other  // Data type of the std::array
        , size_t Dim >    // Dimension of the std::array
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( const std::array<Other,Dim>& array )
{
   resize( Dim, false );

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DynamicVector.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
//
// The vector is resized according to the given N-dimensional vector and initialized as a
// copy of this vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( const DynamicVector& rhs )
{
   if( &rhs == this ) return *this;

   resize( rhs.size_, false );
   smpAssign( *this, ~rhs );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for DynamicVector.
//
// \param rhs The vector to be moved into this instance.
// \return Reference to the assigned vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( DynamicVector&& rhs ) noexcept
{
   deallocate( v_ );

   size_     = rhs.size_;
   capacity_ = rhs.capacity_;
   v_        = rhs.v_;

   rhs.size_     = 0UL;
   rhs.capacity_ = 0UL;
   rhs.v_        = nullptr;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
//
// The vector is resized according to the given vector and initialized as a copy of this vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator=( const Vector<VT,TF>& rhs )
{
   if( (~rhs).canAlias( this ) ) {
      DynamicVector tmp( ~rhs );
      swap( tmp );
   }
   else {
      resize( (~rhs).size(), false );
      if( IsSparseVector_v<VT> )
         reset();
      smpAssign( *this, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator+=( const Vector<VT,TF>& rhs )
{
   if( (~rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (~rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( ~rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator-=( const Vector<VT,TF>& rhs )
{
   if( (~rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (~rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( ~rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator*=( const Vector<VT,TF>& rhs )
{
   if( (~rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( IsSparseVector_v<VT> || (~rhs).canAlias( this ) ) {
      DynamicVector<Type,TF> tmp( *this * (~rhs) );
      swap( tmp );
   }
   else {
      smpMultAssign( *this, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator/=( const DenseVector<VT,TF>& rhs )
{
   if( (~rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (~rhs).canAlias( this ) ) {
      DynamicVector<Type,TF> tmp( *this / (~rhs) );
      swap( tmp );
   }
   else {
      smpDivAssign( *this, ~rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the vector.
// \exception std::invalid_argument Invalid vector size for cross product.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::operator%=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< This, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size_ != 3UL || (~rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (~rhs) );
   assign( *this, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the vector.
//
// \return The size of the vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t DynamicVector<Type,TF>::size() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the minimum capacity of the vector.
//
// \return The minimum capacity of the vector.
//
// This function returns the minimum capacity of the vector, which corresponds to the current
// size plus padding.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t DynamicVector<Type,TF>::spacing() const noexcept
{
   return addPadding( size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the vector.
//
// \return The maximum capacity of the vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t DynamicVector<Type,TF>::capacity() const noexcept
{
   return capacity_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the vector.
//
// \return The number of non-zero elements in the vector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t DynamicVector<Type,TF>::nonZeros() const
{
   size_t nonzeros( 0 );

   for( size_t i=0UL; i<size_; ++i ) {
      if( !isDefault( v_[i] ) )
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::reset()
{
   using blaze::clear;
   for( size_t i=0UL; i<size_; ++i )
      clear( v_[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the vector.
//
// \return void
//
// After the clear() function, the size of the vector is 0.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::clear()
{
   resize( 0UL, false );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
//
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function resizes the vector using the given size to \a n. During this operation, new
// dynamic memory may be allocated in case the capacity of the vector is too small. Note that
// this function may invalidate all existing views (subvectors, ...) on the vector if it is
// used to shrink the vector. Additionally, the resize operation potentially changes all vector
// elements. In order to preserve the old vector values, the \a preserve flag can be set to
// \a true. However, new vector elements are not initialized!
//
// The following example illustrates the resize operation of a vector of size 2 to a vector of
// size 4. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::resize( size_t n, bool preserve )
{
   using std::swap;

   if( n > capacity_ )
   {
      // Allocating a new array
      const size_t newCapacity( addPadding( n ) );
      Type* BLAZE_RESTRICT tmp = allocate<Type>( newCapacity );

      // Initializing the new array
      if( preserve ) {
         transfer( v_, v_+size_, tmp );
      }

      if( IsVectorizable_v<Type> ) {
         for( size_t i=size_; i<newCapacity; ++i )
            tmp[i] = Type();
      }

      // Replacing the old array
      swap( v_, tmp );
      deallocate( tmp );
      capacity_ = newCapacity;
   }
   else if( IsVectorizable_v<Type> && n < size_ )
   {
      for( size_t i=n; i<size_; ++i )
         v_[i] = Type();
   }

   size_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the vector.
//
// \param n Number of additional vector elements.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function increases the vector size by \a n elements. During this operation, new dynamic
// memory may be allocated in case the capacity of the vector is too small. Therefore this
// function potentially changes all vector elements. In order to preserve the old vector values,
// the \a preserve flag can be set to \a true. However, new vector elements are not initialized!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::extend( size_t n, bool preserve )
{
   resize( size_+n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the vector.
//
// \param n The new minimum capacity of the vector.
// \return void
//
// This function increases the capacity of the vector to at least \a n elements. The current
// values of the vector elements are preserved.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::reserve( size_t n )
{
   using std::swap;

   if( n > capacity_ )
   {
      // Allocating a new array
      const size_t newCapacity( addPadding( n ) );
      Type* BLAZE_RESTRICT tmp = allocate<Type>( newCapacity );

      // Initializing the new array
      transfer( v_, v_+size_, tmp );

      if( IsVectorizable_v<Type> ) {
         for( size_t i=size_; i<newCapacity; ++i )
            tmp[i] = Type();
      }

      // Replacing the old array
      swap( tmp, v_ );
      deallocate( tmp );
      capacity_ = newCapacity;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the vector by removing unused capacity. Please note
// that due to padding the capacity might not be reduced exactly to size(). Please also note
// that in case a reallocation occurs, all iterators (including end() iterators), all pointers
// and references to elements of this vector are invalidated.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::shrinkToFit()
{
   if( spacing() < capacity_ ) {
      DynamicVector( *this ).swap( *this );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
//
// \param v The vector to be swapped.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void DynamicVector<Type,TF>::swap( DynamicVector& v ) noexcept
{
   using std::swap;

   swap( size_, v.size_ );
   swap( capacity_, v.capacity_ );
   swap( v_, v.v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Add the necessary amount of padding to the given value.
//
// \param value The value to be padded.
// \return The padded value.
//
// This function increments the given \a value by the necessary amount of padding based on the
// vector's data type \a Type.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t DynamicVector<Type,TF>::addPadding( size_t value ) const noexcept
{
   if( IsVectorizable_v<Type> )
      return nextMultiple<size_t>( value, SIMDSIZE );
   else return value;
}
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scaling of the vector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the vector scaling.
// \return Reference to the vector.
//
// This function scales the vector by applying the given scalar value \a scalar to each element
// of the vector. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   blaze::DynamicVector<int> a;
   // ... Initialization
   a *= 4;        // Scaling of the vector
   a.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the scalar value
inline DynamicVector<Type,TF>& DynamicVector<Type,TF>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] *= scalar;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  DEBUGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the invariants of the dynamic vector are intact.
//
// \return \a true in case the dynamic vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dynamic vector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline bool DynamicVector<Type,TF>::isIntact() const noexcept
{
   if( size_ > capacity_ )
      return false;

   if( IsVectorizable_v<Type> ) {
      for( size_t i=size_; i<capacity_; ++i ) {
         if( v_[i] != Type() )
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
/*!\brief Returns whether the vector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool DynamicVector<Type,TF>::canAlias( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address is aliased with the vector. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool DynamicVector<Type,TF>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector is properly aligned in memory.
//
// \return \a true in case the vector is aligned, \a false if not.
//
// This function returns whether the vector is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of the vector are guaranteed to conform to the alignment
// restrictions of the element type \a Type.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline bool DynamicVector<Type,TF>::isAligned() const noexcept
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector can be used in SMP assignments.
//
// \return \a true in case the vector can be used in SMP assignments, \a false if not.
//
// This function returns whether the vector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// vector).
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline bool DynamicVector<Type,TF>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense vector. The index
// must be smaller than the number of vector elements and it must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE typename DynamicVector<Type,TF>::SIMDType
   DynamicVector<Type,TF>::load( size_t index ) const noexcept
{
   return loada( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE typename DynamicVector<Type,TF>::SIMDType
   DynamicVector<Type,TF>::loada( size_t index ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid alignment detected" );

   return loada( v_+index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE typename DynamicVector<Type,TF>::SIMDType
   DynamicVector<Type,TF>::loadu( size_t index ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );

   return loadu( v_+index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense vector. The index
// must be smaller than the number of vector elements and it must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE void
   DynamicVector<Type,TF>::store( size_t index, const SIMDType& value ) noexcept
{
   storea( index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE void
   DynamicVector<Type,TF>::storea( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid alignment detected" );

   storea( v_+index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense vector.
// The index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE void
   DynamicVector<Type,TF>::storeu( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );

   storeu( v_+index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense vector. The index must be smaller than the number of vector elements and it must be
// a multiple of the number of values inside the SIMD element. This function must \b NOT be
// called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
BLAZE_ALWAYS_INLINE void
   DynamicVector<Type,TF>::stream( size_t index, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid alignment detected" );

   stream( v_+index, value );
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::assign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] = (~rhs)[i    ];
      v_[i+1UL] = (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      v_[ipos] = (~rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::assign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i=0UL;
   Iterator left( begin() );
   ConstIterator_t<VT> right( (~rhs).begin() );

   if( useStreaming &&
       ( size_ > ( cacheSize/( sizeof(Type) * 3UL ) ) ) &&
       !(~rhs).isAliased( this ) )
   {
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<size_; ++i ) {
         *left = *right; ++left; ++right;
      }
   }
   else
   {
      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<size_; ++i ) {
         *left = *right; ++left; ++right;
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DynamicVector<Type,TF>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   for( auto element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      v_[element->index()] = element->value();
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::addAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] += (~rhs)[i    ];
      v_[i+1UL] += (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      v_[ipos] += (~rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::addAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (~rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && i<size_; ++i ) {
      *left += *right; ++left; ++right;
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DynamicVector<Type,TF>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   for( auto element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      v_[element->index()] += element->value();
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::subAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] -= (~rhs)[i    ];
      v_[i+1UL] -= (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      v_[ipos] -= (~rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::subAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (~rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && i<size_; ++i ) {
      *left -= *right; ++left; ++right;
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DynamicVector<Type,TF>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   for( auto element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      v_[element->index()] -= element->value();
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::multAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] *= (~rhs)[i    ];
      v_[i+1UL] *= (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      v_[ipos] *= (~rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::multAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (~rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; remainder && i<size_; ++i ) {
      *left *= *right; ++left; ++right;
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
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void DynamicVector<Type,TF>::multAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const DynamicVector tmp( serial( *this ) );

   reset();

   for( auto element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      v_[element->index()] = tmp[element->index()] * element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisior.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::divAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] /= (~rhs)[i    ];
      v_[i+1UL] /= (~rhs)[i+1UL];
   }
   if( ipos < (~rhs).size() )
      v_[ipos] /= (~rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline auto DynamicVector<Type,TF>::divAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT> right( (~rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size_; ++i ) {
      *left /= *right; ++left; ++right;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DYNAMICVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DynamicVector operators */
//@{
template< typename Type, bool TF >
void reset( DynamicVector<Type,TF>& v );

template< typename Type, bool TF >
void clear( DynamicVector<Type,TF>& v );

template< RelaxationFlag RF, typename Type, bool TF >
bool isDefault( const DynamicVector<Type,TF>& v );

template< typename Type, bool TF >
bool isIntact( const DynamicVector<Type,TF>& v ) noexcept;

template< typename Type, bool TF >
void swap( DynamicVector<Type,TF>& a, DynamicVector<Type,TF>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dynamic vector.
// \ingroup dynamic_vector
//
// \param v The dynamic vector to be resetted.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void reset( DynamicVector<Type,TF>& v )
{
   v.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dynamic vector.
// \ingroup dynamic_vector
//
// \param v The dynamic vector to be cleared.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void clear( DynamicVector<Type,TF>& v )
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dynamic vector is in default state.
// \ingroup dynamic_vector
//
// \param v The dynamic vector to be tested for its default state.
// \return \a true in case the given vector's size is zero, \a false otherwise.
//
// This function checks whether the dynamic vector is in default (constructed) state, i.e. if
// it's size is 0. In case it is in default state, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isDefault() function:

   \code
   blaze::DynamicVector<int> a;
   // ... Resizing and initialization
   if( isDefault( a ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( a ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type      // Data type of the vector
        , bool TF >          // Transpose flag
inline bool isDefault( const DynamicVector<Type,TF>& v )
{
   return ( v.size() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given dynamic vector are intact.
// \ingroup dynamic_vector
//
// \param v The dynamic vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dynamic vector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::DynamicVector<int> a;
   // ... Resizing and initialization
   if( isIntact( a ) ) { ... }
   \endcode
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline bool isIntact( const DynamicVector<Type,TF>& v ) noexcept
{
   return v.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
// \ingroup dynamic_vector
//
// \param a The first vector to be swapped.
// \param b The second vector to be swapped.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void swap( DynamicVector<Type,TF>& a, DynamicVector<Type,TF>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF >
struct HasConstDataAccess< DynamicVector<T,TF> >
   : public TrueType
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
template< typename T, bool TF >
struct HasMutableDataAccess< DynamicVector<T,TF> >
   : public TrueType
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
template< typename T, bool TF >
struct IsAligned< DynamicVector<T,TF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGUOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF >
struct IsContiguous< DynamicVector<T,TF> >
   : public TrueType
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
template< typename T, bool TF >
struct IsPadded< DynamicVector<T,TF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF >
struct IsResizable< DynamicVector<T,TF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSHRINKABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF >
struct IsShrinkable< DynamicVector<T,TF> >
   : public TrueType
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
template< typename T1, typename T2 >
struct AddTraitEval2< T1, T2
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  ( IsDenseVector_v<T1> || IsDenseVector_v<T2> ) &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< AddTrait_t<ET1,ET2>, TransposeFlag_v<T1> >;
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
template< typename T1, typename T2 >
struct SubTraitEval2< T1, T2
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  ( IsDenseVector_v<T1> || IsDenseVector_v<T2> ) &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< SubTrait_t<ET1,ET2>, TransposeFlag_v<T1> >;
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
template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsDenseVector_v<T1> &&
                                   IsNumeric_v<T2> &&
                                   ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) > >
{
   using ET1 = ElementType_t<T1>;

   using Type = DynamicVector< MultTrait_t<ET1,T2>, TransposeFlag_v<T1> >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsNumeric_v<T1> &&
                                   IsDenseVector_v<T2> &&
                                   ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) > >
{
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< MultTrait_t<T1,ET2>, TransposeFlag_v<T2> >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                     ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                   IsDenseVector_v<T1> &&
                                   IsDenseVector_v<T2> &&
                                   ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                   ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                   ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< MultTrait_t<ET1,ET2>, TransposeFlag_v<T1> >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsColumnVector_v<T2> &&
                                   ( IsDenseMatrix_v<T1> || IsDenseVector_v<T2> ) &&
                                   ( Size_v<T1,0UL> == DefaultSize_v &&
                                     ( !IsSquare_v<T1> || Size_v<T2,0UL> == DefaultSize_v ) ) &&
                                   ( MaxSize_v<T1,0UL> == DefaultMaxSize_v &&
                                     ( !IsSquare_v<T1> || MaxSize_v<T2,0UL> == DefaultMaxSize_v ) ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< MultTrait_t<ET1,ET2>, false >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsRowVector_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsDenseVector_v<T1> || IsDenseMatrix_v<T2> ) &&
                                   ( Size_v<T2,1UL> == DefaultSize_v &&
                                     ( !IsSquare_v<T2> || Size_v<T1,0UL> == DefaultSize_v ) ) &&
                                   ( MaxSize_v<T2,1UL> == DefaultMaxSize_v &&
                                     ( !IsSquare_v<T2> || MaxSize_v<T1,0UL> == DefaultMaxSize_v ) ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< MultTrait_t<ET1,ET2>, true >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  KRONTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct KronTraitEval2< T1, T2
                     , EnableIf_t< IsDenseVector_v<T1> &&
                                   IsDenseVector_v<T2> &&
                                   ( ( Size_v<T1,0UL> == DefaultSize_v ) ||
                                     ( Size_v<T2,0UL> == DefaultSize_v ) ) &&
                                   ( ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) ||
                                     ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< MultTrait_t<ET1,ET2>, TransposeFlag_v<T2> >;
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
template< typename T1, typename T2 >
struct DivTraitEval2< T1, T2
                    , EnableIf_t< IsDenseVector_v<T1> &&
                                  IsNumeric_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) > >
{
   using ET1 = ElementType_t<T1>;

   using Type = DynamicVector< DivTrait_t<ET1,T2>, TransposeFlag_v<T1> >;
};

template< typename T1, typename T2 >
struct DivTraitEval2< T1, T2
                    , EnableIf_t< IsDenseVector_v<T1> &&
                                  IsDenseVector_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< DivTrait_t<ET1,ET2>, TransposeFlag_v<T1> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP >
struct UnaryMapTraitEval2< T, OP
                         , EnableIf_t< IsDenseVector_v<T> &&
                                       Size_v<T,0UL> == DefaultSize_v &&
                                       MaxSize_v<T,0UL> == DefaultMaxSize_v > >
{
   using ET = ElementType_t<T>;

   using Type = DynamicVector< MapTrait_t<ET,OP>, TransposeFlag_v<T> >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval2< T1, T2, OP
                          , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                          ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                        Size_v<T1,0UL> == DefaultSize_v &&
                                        Size_v<T2,0UL> == DefaultSize_v &&
                                        MaxSize_v<T1,0UL> == DefaultMaxSize_v &&
                                        MaxSize_v<T2,0UL> == DefaultMaxSize_v > >
{
   using ET1 = ElementType_t<T1>;
   using ET2 = ElementType_t<T2>;

   using Type = DynamicVector< MapTrait_t<ET1,ET2,OP>, TransposeFlag_v<T1> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REDUCETRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP, ReductionFlag RF >
struct PartialReduceTraitEval2< T, OP, RF
                              , EnableIf_t< IsMatrix_v<T> &&
                                            ( Size_v<T,0UL> == DefaultSize_v ||
                                              Size_v<T,1UL> == DefaultSize_v ) &&
                                            ( MaxSize_v<T,0UL> == DefaultMaxSize_v ||
                                              MaxSize_v<T,1UL> == DefaultMaxSize_v ) > >
{
   using ET = ElementType_t<T>;
   using RT = decltype( std::declval<OP>()( std::declval<ET>(), std::declval<ET>() ) );

   static constexpr bool TF = ( RF == columnwise );

   using Type = DynamicVector<RT,TF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SOLVETRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SolveTraitEval2< T1, T2
                      , EnableIf_t< IsDenseMatrix_v<T1> &&
                                    IsDenseVector_v<T2> &&
                                    ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                    ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                    ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                    ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) > >
{
   using Type = DynamicVector< ElementType_t<T2>, TransposeFlag_v<T2> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool TF, typename T2 >
struct HighType< DynamicVector<T1,TF>, DynamicVector<T2,TF> >
{
   using Type = DynamicVector< typename HighType<T1,T2>::Type, TF >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool TF, typename T2 >
struct LowType< DynamicVector<T1,TF>, DynamicVector<T2,TF> >
{
   using Type = DynamicVector< typename LowType<T1,T2>::Type, TF >;
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
template< typename VT >
struct SubvectorTraitEval2< VT, inf, inf
                          , EnableIf_t< IsDenseVector_v<VT> &&
                                        Size_v<VT,0UL> == DefaultSize_v &&
                                        MaxSize_v<VT,0UL> == DefaultMaxSize_v > >
{
   using Type = DynamicVector< RemoveConst_t< ElementType_t<VT> >, TransposeFlag_v<VT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ELEMENTSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT >
struct ElementsTraitEval2< VT, 0UL
                         , EnableIf_t< IsDenseVector_v<VT> &&
                                       Size_v<VT,0UL> == DefaultSize_v &&
                                       MaxSize_v<VT,0UL> == DefaultMaxSize_v > >
{
   using Type = DynamicVector< RemoveConst_t< ElementType_t<VT> >, TransposeFlag_v<VT> >;
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
template< typename MT, size_t I >
struct RowTraitEval2< MT, I
                    , EnableIf_t< IsDenseMatrix_v<MT> &&
                                  Size_v<MT,1UL> == DefaultSize_v &&
                                  MaxSize_v<MT,1UL> == DefaultMaxSize_v > >
{
   using Type = DynamicVector< RemoveConst_t< ElementType_t<MT> >, true >;
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
template< typename MT, size_t I >
struct ColumnTraitEval2< MT, I
                       , EnableIf_t< IsDenseMatrix_v<MT> &&
                                     Size_v<MT,0UL> == DefaultSize_v &&
                                     MaxSize_v<MT,0UL> == DefaultMaxSize_v > >
{
   using Type = DynamicVector< RemoveConst_t< ElementType_t<MT> >, false >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, ptrdiff_t I >
struct BandTraitEval2< MT, I
                     , EnableIf_t< IsDenseMatrix_v<MT> &&
                                   ( Size_v<MT,0UL> == DefaultSize_v ||
                                     Size_v<MT,1UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<MT,0UL> == DefaultMaxSize_v ||
                                     MaxSize_v<MT,1UL> == DefaultMaxSize_v ) > >
{
   using Type = DynamicVector< RemoveConst_t< ElementType_t<MT> >, defaultTransposeFlag >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
