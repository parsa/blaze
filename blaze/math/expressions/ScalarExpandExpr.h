//=================================================================================================
/*!
//  \file blaze/math/expressions/ScalarExpandExpr.h
//  \brief Header file for the scalar expansion expression
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SCALAREXPANDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SCALAREXPANDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/DVecExpandExpr.h>
#include <blaze/math/expressions/ExpandExpr.h>
#include <blaze/math/expressions/ExpandExprData.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/Transformation.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/ExpandTrait.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/views/elements/ElementsData.h>
#include <blaze/math/views/subvector/SubvectorData.h>
#include <blaze/system/Inline.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SCALAREXPANDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for scalar expansion.
// \ingroup scalar_expression
//
// The ScalarExpandExpr class represents the compile time expression for expansions of scalars.
*/
template< typename ST       // Type of the scalar
        , bool TF           // Transpose flag
        , size_t... CEAs >  // Compile time expansion arguments
class ScalarExpandExpr
   : public ExpandExpr< DenseVector< ScalarExpandExpr<ST,TF,CEAs...>, TF > >
   , private Transformation
   , private ExpandExprData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ExpandExprData<CEAs...>;  //!< The type of the ExpandExprData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This          = ScalarExpandExpr<ST,TF,CEAs...>;  //!< Type of this DVecExpandExpr instance.
   using BaseType      = DenseVector<This,TF>;             //!< Base type of this ScalarExpandExpr instance.
   using ResultType    = ExpandTrait_t<ST,CEAs...>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;      //!< Transpose type for expression template evaluations.
   using ElementType   = ST;                               //!< Resulting element type.
   using ReturnType    = const ST&;                        //!< Return type for expression template evaluations.
   using CompositeType = const ScalarExpandExpr&;          //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::random_access_iterator_tag;  //!< The iterator category.
      using ValueType        = ElementType;                      //!< Type of the underlying elements.
      using PointerType      = ElementType*;                     //!< Pointer return type.
      using ReferenceType    = ElementType&;                     //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                        //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param index Index to the initial element.
      // \param scalar Scalar of the expansion expression.
      */
      explicit inline ConstIterator( size_t index, ST scalar )
         : index_ ( index  )  // Index to the current element
         , scalar_( scalar )  // Scalar of the expansion expression
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ConstIterator& operator+=( size_t inc ) {
         index_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline ConstIterator& operator-=( size_t dec ) {
         index_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ConstIterator& operator++() {
         ++index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator++( int ) {
         return ConstIterator( index_++, scalar_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ConstIterator& operator--() {
         --index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator--( int ) {
         return ConstIterator( index_--, scalar_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReturnType operator*() const {
         return scalar_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Access to the SIMD elements of the vector.
      //
      // \return The resulting SIMD element.
      */
      inline auto load() const noexcept {
         return set( scalar_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ConstIterator& rhs ) const {
         return index_ < rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ConstIterator& rhs ) const {
         return index_ > rhs.index_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ConstIterator& rhs ) const {
         return index_ <= rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ConstIterator& rhs ) const {
         return index_ >= rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return index_ - rhs.index_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ConstIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( const ConstIterator& it, size_t inc ) {
         return ConstIterator( it.index_ + inc, it.scalar_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ConstIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( size_t inc, const ConstIterator& it ) {
         return ConstIterator( it.index_ + inc, it.scalar_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ConstIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ConstIterator operator-( const ConstIterator& it, size_t dec ) {
         return ConstIterator( it.index_ - dec, it.scalar_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      size_t index_;   //!< Index to the current element.
      ST     scalar_;  //!< Scalar of the expansion expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = true;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = true;
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the ScalarExpandExpr class.
   //
   // \param scalar The scalar value of the expansion expression.
   // \param args The runtime expansion expression arguments.
   */
   template< typename... REAs >  // Runtime expansion arguments
   explicit inline ScalarExpandExpr( ST scalar, REAs... args ) noexcept
      : DataType( args... )  // Base class initialization
      , scalar_ ( scalar )   // Scalar value of the expansion expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      MAYBE_UNUSED( index );
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return scalar_;
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

   //**Load function*******************************************************************************
   /*!\brief Access to the SIMD elements of the vector.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE auto load( size_t index ) const noexcept {
      MAYBE_UNUSED( index );
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
      return set( scalar_ );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( 0UL, scalar_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( size(), scalar_ );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return expansion();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the scalar value.
   //
   // \return The scalar value.
   */
   inline ST operand() const noexcept {
      return scalar_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::expansion;
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      MAYBE_UNUSED( alias );
      return false;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      MAYBE_UNUSED( alias );
      return false;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return true;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return true;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   ST scalar_;  //!< Scalar value of the expansion expression.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param expansion The expansion.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::columnVector;

   int scalar = 5;

   blaze::DynamicVector<int,columnVector> v;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3-dimensional column vector
   //
   //    ( 5 )
   //    ( 5 )
   //    ( 5 )
   //
   v = expand( scalar, 3UL );
   \endcode
*/
template< typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expand( ST scalar, size_t expansion )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const ScalarExpandExpr<ST,defaultTransposeFlag>;
   return ReturnType( scalar, expansion );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::columnVector;

   int scalar = 5;

   blaze::DynamicVector<int,columnVector> v;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3-dimensional column vector
   //
   //    ( 5 )
   //    ( 5 )
   //    ( 5 )
   //
   v = expand<3UL>( scalar );
   \endcode
*/
template< size_t E     // Compile time expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expand( ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const ScalarExpandExpr<ST,defaultTransposeFlag,E>;
   return ReturnType( scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param expansion The expansion.
// \return The expansion of the scalar value.
//
// This auxiliary overload of the \c expand() function accepts both a compile time and a runtime
// expansion. The runtime argument is discarded in favor of the compile time argument.
*/
template< size_t E     // Compile time expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expand( ST scalar, size_t expansion )
{
   MAYBE_UNUSED( expansion );

   BLAZE_FUNCTION_TRACE;

   using ReturnType = const ScalarExpandExpr<ST,defaultTransposeFlag,E>;
   return ReturnType( scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param expansion The expansion.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::rowVector;

   int scalar = 5;

   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3-dimensional row vector
   //
   //    ( 5  5  5 )
   //
   v = expandTo<rowVector>( scalar, 3UL );
   \endcode
*/
template< bool TF      // Target transpose flag
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expandTo( ST scalar, size_t expansion )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const ScalarExpandExpr<ST,TF>;
   return ReturnType( scalar, expansion );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::rowVector;

   int scalar = 5;

   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3-dimensional row vector
   //
   //    ( 5  5  5 )
   //
   v = expandTo<rowVector,3UL>( scalar );
   \endcode
*/
template< bool TF      // Target transpose flag
        , size_t E     // Compile time expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expandTo( ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const ScalarExpandExpr<ST,TF,E>;
   return ReturnType( scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param expansion The expansion.
// \return The expansion of the scalar value.
//
// This auxiliary overload of the \c expandTo() function accepts both a compile time and a runtime
// expansion. The runtime argument is discarded in favor of the compile time argument.
*/
template< bool TF      // Target transpose flag
        , size_t E     // Compile time expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expandTo( ST scalar, size_t expansion )
{
   MAYBE_UNUSED( expansion );

   BLAZE_FUNCTION_TRACE;

   using ReturnType = const ScalarExpandExpr<ST,TF,E>;
   return ReturnType( scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param rowExpansion The row expansion.
// \param columnExpansion The column expansion.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::rowMajor;

   int scalar = 5;

   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3x5 row-major matrix
   //
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //
   A = expand( scalar, 3UL, 5UL );
   \endcode
*/
template< typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expand( ST scalar, size_t rowExpansion, size_t columnExpansion )
{
   BLAZE_FUNCTION_TRACE;

   const size_t m( defaultStorageOrder ? columnExpansion : rowExpansion );
   const size_t n( defaultStorageOrder ? rowExpansion : columnExpansion );

   return expand( expandTo<!defaultStorageOrder>( scalar, n ), m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::rowMajor;

   int scalar = 5;

   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3x5 row-major matrix
   //
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //
   A = expand<3UL,5UL>( scalar );
   \endcode
*/
template< size_t R     // Compile time row expansion argument
        , size_t C     // Compile time column expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expand( ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   constexpr size_t M( defaultStorageOrder ? C : R );
   constexpr size_t N( defaultStorageOrder ? R : C );

   return expand<M>( expandTo<!defaultStorageOrder,N>( scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param rowExpansion The row expansion.
// \param rowExpansion The column expansion.
// \return The expansion of the scalar value.
//
// This auxiliary overload of the \c expand() function accepts both a compile time and a runtime
// expansions. The runtime arguments are discarded in favor of the compile time arguments.
*/
template< size_t R     // Compile time row expansion argument
        , size_t C     // Compile time column expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expand( ST scalar, size_t rowExpansion, size_t columnExpansion )
{
   MAYBE_UNUSED( rowExpansion, columnExpansion );

   BLAZE_FUNCTION_TRACE;

   constexpr size_t M( defaultStorageOrder ? C : R );
   constexpr size_t N( defaultStorageOrder ? R : C );

   return expand<M>( expandTo<!defaultStorageOrder,N>( scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param rowExpansion The row expansion.
// \param columnExpansion The column expansion.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::columnMajor;

   int scalar = 5;

   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3x5 column-major matrix
   //
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //
   A = expandTo<columnMajor>( scalar, 3UL, 5UL );
   \endcode
*/
template< bool SO      // Target storage order
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expandTo( ST scalar, size_t rowExpansion, size_t columnExpansion )
{
   BLAZE_FUNCTION_TRACE;

   const size_t m( SO ? columnExpansion : rowExpansion );
   const size_t n( SO ? rowExpansion : columnExpansion );

   return expand( expandTo<!SO>( scalar, n ), m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \return The expansion of the scalar value.
//
// This function returns an expression representing the expansion of the given scalar value:

   \code
   using blaze::columnMajor;

   int scalar = 5;

   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization

   // Expansion of 'scalar' to a 3x5 column-major matrix
   //
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //    ( 5  5  5  5  5 )
   //
   A = expandTo<columnMajor,3UL,5UL>( scalar );
   \endcode
*/
template< bool SO      // Target storage order
        , size_t R     // Compile time row expansion argument
        , size_t C     // Compile time column expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expandTo( ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   constexpr size_t M( SO ? C : R );
   constexpr size_t N( SO ? R : C );

   return expand<M>( expandTo<!SO,N>( scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expansion of the given scalar value.
// \ingroup scalar
//
// \param scalar The scalar value to be expanded.
// \param rowExpansion The row expansion.
// \param rowExpansion The column expansion.
// \return The expansion of the scalar value.
//
// This auxiliary overload of the \c expandTo() function accepts both a compile time and a runtime
// expansions. The runtime arguments are discarded in favor of the compile time arguments.
*/
template< bool SO      // Target storage order
        , size_t R     // Compile time row expansion argument
        , size_t C     // Compile time column expansion argument
        , typename ST  // Type of the scalar
        , EnableIf_t< IsNumeric_v<ST> >* = nullptr >
inline decltype(auto) expandTo( ST scalar, size_t rowExpansion, size_t columnExpansion )
{
   MAYBE_UNUSED( rowExpansion, columnExpansion );

   BLAZE_FUNCTION_TRACE;

   constexpr size_t M( SO ? C : R );
   constexpr size_t N( SO ? R : C );

   return expand<M>( expandTo<!SO,N>( scalar ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (SUBVECTOR)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given scalar expansion operation.
// \ingroup subvector
//
// \param vector The constant scalar expansion operation.
// \param args The runtime subvector arguments.
// \return View on the specified subvector of the expansion operation.
//
// This function returns an expression representing the specified subvector of the given scalar
// expansion operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time subvector arguments
        , typename ST         // Type of the scalar
        , bool TF             // Transpose flag
        , size_t... CEAs      // Compile time expansion arguments
        , typename... RSAs >  // Runtime subvector arguments
inline decltype(auto) subvector( const ScalarExpandExpr<ST,TF,CEAs...>& vector, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const SubvectorData<CSAs...> sd( args... );

   return expandTo<TF>( vector.operand(), sd.size() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ELEMENTS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given scalar expansion operation.
// \ingroup elements
//
// \param vector The constant scalar expansion operation.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the expansion operation.
//
// This function returns an expression representing the specified selection of elements on the
// given scalar expansion operation.
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename ST         // Type of the scalar
        , bool TF             // Transpose flag
        , size_t... CEAs      // Compile time expansion arguments
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const ScalarExpandExpr<ST,TF,CEAs...>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   return expandTo<TF>( vector.operand(), sizeof...(Is)+1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given scalar expansion operation.
// \ingroup elements
//
// \param vector The constant scalar expansion operation.
// \param indices Pointer to the first index of the selected elements.
// \param n The total number of indices.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the expansion operation.
//
// This function returns an expression representing the specified selection of elements on the
// given scalar expansion operation.
*/
template< typename ST         // Type of the scalar
        , bool TF             // Transpose flag
        , size_t... CEAs      // Compile time expansion arguments
        , typename T          // Type of the element indices
        , typename... REAs >  // Runtime element arguments
inline decltype(auto)
   elements( const ScalarExpandExpr<ST,TF,CEAs...>& vector, T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( indices, args... );

   return expandTo<TF>( vector.operand(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given scalar expansion operation.
// \ingroup elements
//
// \param vector The constant scalar expansion operation.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the expansion operation.
//
// This function returns an expression representing the specified selection of elements on the
// given scalar expansion operation.
*/
template< typename ST         // Type of the scalar
        , bool TF             // Transpose flag
        , size_t... CEAs      // Compile time expansion arguments
        , typename P          // Type of the index producer
        , typename... REAs >  // Runtime element arguments
inline decltype(auto)
   elements( const ScalarExpandExpr<ST,TF,CEAs...>& vector, P p, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( p, args... );

   return expandTo<TF>( vector.operand(), n );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename ST, bool TF, size_t... CEAs >
struct IsAligned< ScalarExpandExpr<ST,TF,CEAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
