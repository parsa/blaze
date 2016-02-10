//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecScalarDivExpr.h
//  \brief Header file for the dense vector/scalar division expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSCALARDIVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSCALARDIVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecScalarDivExpr.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsMultExpr.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSCALARDIVEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for divisions of a dense vector by a scalar.
// \ingroup dense_vector_expression
//
// The DVecScalarDivExpr class represents the compile time expression for divisions of dense
// vectors by scalar values.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar value
        , bool TF >    // Transpose flag
class DVecScalarDivExpr : public DenseVector< DVecScalarDivExpr<VT,ST,TF>, TF >
                        , private VecScalarDivExpr
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     RT;  //!< Result type of the dense vector expression.
   typedef typename VT::ReturnType     RN;  //!< Return type of the dense vector expression.
   typedef typename VT::ElementType    ET;  //!< Element type of the dense vector expression.
   typedef typename VT::CompositeType  CT;  //!< Composite type of the dense vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If the vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   enum { returnExpr = !IsTemporary<RN>::value };

   //! Expression return type for the subscript operator.
   typedef typename DivExprTrait<RN,ST>::Type  ExprReturnType;
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the division expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the division expression. In case the given dense
       vector expression of type \a VT is a computation expression and requires an intermediate
       evaluation, \a useAssign will be set to 1 and the division expression will be evaluated
       via the \a assign function family. Otherwise \a useAssign will be set to 0 and the
       expression will be evaluated via the subscript operator. */
   enum { useAssign = IsComputation<VT>::value && RequiresEvaluation<VT>::value };

   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct UseAssign {
      enum { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The UseSMPAssign struct is a helper struct for the selection of the parallel evaluation
       strategy. In case either the target vector or the dense vector operand is not SMP assignable
       and the vector operand is a computation expression that requires an intermediate evaluation,
       \a value is set to 1 and the expression specific evaluation strategy is selected. Otherwise
       \a value is set to 0 and the default strategy is chosen. */
   template< typename VT2 >
   struct UseSMPAssign {
      enum { value = ( !VT2::smpAssignable || !VT::smpAssignable ) && useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DVecScalarDivExpr<VT,ST,TF>                 This;           //!< Type of this DVecScalarDivExpr instance.
   typedef typename DivTrait<RT,ST>::Type              ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DVecScalarDivExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  LeftOperand;

   //! Composite type of the right-hand side scalar value.
   typedef ST  RightOperand;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::random_access_iterator_tag  IteratorCategory;  //!< The iterator category.
      typedef ElementType                      ValueType;         //!< Type of the underlying elements.
      typedef ElementType*                     PointerType;       //!< Pointer return type.
      typedef ElementType&                     ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                        DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.

      //! ConstIterator type of the dense vector expression.
      typedef typename VT::ConstIterator  IteratorType;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param scalar Scalar of the division expression.
      */
      explicit inline ConstIterator( IteratorType iterator, RightOperand scalar )
         : iterator_( iterator )  // Iterator to the current element
         , scalar_  ( scalar   )  // Scalar of the division expression
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ConstIterator& operator+=( size_t inc ) {
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
      inline ConstIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ConstIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator++( int ) {
         return ConstIterator( iterator_++ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ConstIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator--( int ) {
         return ConstIterator( iterator_-- );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReturnType operator*() const {
         return *iterator_ / scalar_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Access to the intrinsic elements of the vector.
      //
      // \return The resulting intrinsic value.
      */
      inline IntrinsicType load() const {
         return iterator_.load() / set( scalar_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ConstIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ConstIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ConstIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ConstIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
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
         return ConstIterator( it.iterator_ + inc );
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
         return ConstIterator( it.iterator_ + inc );
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
         return ConstIterator( it.iterator_ - dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;  //!< Iterator to the current element.
      RightOperand scalar_;    //!< Scalar of the division expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT::vectorizable &&
                         IsNumeric<ET>::value &&
                         ( IsSame<ET,RightOperand>::value ||
                           IsSame<typename UnderlyingElement<ET>::Type,RightOperand>::value ) &&
                         IntrinsicTrait<RightOperand>::division };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = VT::smpAssignable };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecScalarDivExpr class.
   //
   // \param vector The left-hand side dense vector of the division expression.
   // \param scalar The right-hand side scalar of the division expression.
   */
   explicit inline DVecScalarDivExpr( const VT& vector, ST scalar )
      : vector_( vector )  // Left-hand side dense vector of the division expression
      , scalar_( scalar )  // Right-hand side scalar of the division expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < vector_.size(), "Invalid vector access index" );
      return vector_[index] / scalar_;
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
      if( index >= vector_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Load function*******************************************************************************
   /*!\brief Access to the intrinsic elements of the vector.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE IntrinsicType load( size_t index ) const {
      typedef IntrinsicTrait<ElementType>  IT;
      BLAZE_INTERNAL_ASSERT( index < vector_.size() , "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL, "Invalid vector access index" );
      return vector_.load( index ) / set( scalar_ );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( vector_.begin(), scalar_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( vector_.end(), scalar_ );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return vector_.size();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const {
      return vector_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side scalar operand.
   //
   // \return The right-hand side scalar operand.
   */
   inline RightOperand rightOperand() const {
      return scalar_;
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
      return IsComputation<VT>::value && vector_.canAlias( alias );
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

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const {
      return vector_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return vector_.canSMPAssign() || ( size() > SMP_DVECSCALARMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vector_;  //!< Left-hand side dense vector of the division expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the division expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-scalar
   // division expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the vector
   // operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( ~lhs, rhs.vector_ );
      assign( ~lhs, (~lhs) / rhs.scalar_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-scalar division to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-scalar
   // division expression to a sparse vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the vector
   // operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( SparseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( ~lhs, rhs.vector_ );
      (~lhs) /= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // scalar division expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the vector
   // operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      addAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      addAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // scalar division expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the vector operand is
   // a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      subAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      subAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-scalar division expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // vector operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      multAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-
   // scalar division expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT2> >::Type
      smpAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( ~lhs, rhs.vector_ );
      smpAssign( ~lhs, (~lhs) / rhs.scalar_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector-scalar division to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-
   // scalar division expression to a sparse vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseSMPAssign<VT2> >::Type
      smpAssign( SparseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( ~lhs, rhs.vector_ );
      (~lhs) /= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector-scalar division expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT2> >::Type
      smpAddAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAddAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // vector-scalar division expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT2> >::Type
      smpSubAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpSubAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector-scalar division to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side division expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // dense vector-scalar division expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT2> >::Type
      smpMultAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_NOT_BE_FLOATING_POINT_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_NOT_BE_FLOATING_POINT_TYPE( ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Division operator for the divison of a dense vector by a scalar value
//        (\f$ \vec{a}=\vec{b}/s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result vector.
//
// This operator represents the division of a dense vector by a scalar value:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = a / 0.24;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order
// element type of the involved data types \a T1::ElementType and \a T2. Both data types
// \a T1::ElementType and \a T2 have to be supported by the DivTrait class template.
// Note that this operator only works for scalar values of built-in data type.
//
// \note: A division by zero is only checked by an user assert.
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline const typename EnableIf< IsNumeric<T2>, typename DivExprTrait<T1,T2>::Type >::Type
   operator/( const DenseVector<T1,TF>& vec, T2 scalar )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( scalar != T2(0), "Division by zero detected" );

   typedef typename DivExprTrait<T1,T2>::Type  ReturnType;
   typedef typename ReturnType::RightOperand   ScalarType;

   if( IsMultExpr<ReturnType>::value ) {
      return ReturnType( ~vec, ScalarType(1)/ScalarType(scalar) );
   }
   else {
      return ReturnType( ~vec, scalar );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar division
//        expression and a scalar value (\f$ \vec{a}=(\vec{b}/s1)*s2 \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector-scalar division.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// dense vector-scalar division expression and a scalar value.
*/
template< typename VT     // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the dense vector
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< And< IsNumeric<ST2>, IsInvertible<typename DivTrait<ST2,ST1>::Type > >
                              , typename MultExprTrait< DVecScalarDivExpr<VT,ST1,TF>, ST2 >::Type >::Type
   operator*( const DVecScalarDivExpr<VT,ST1,TF>& vec, ST2 scalar )
{
   BLAZE_FUNCTION_TRACE;

   return vec.leftOperand() * ( scalar / vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense vector-
//        scalar division expression (\f$ \vec{a}=s2*(\vec{b}/s1) \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side dense vector-scalar division.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a dense vector-scalar division expression.
*/
template< typename ST1  // Type of the left-hand side scalar
        , typename VT   // Type of the dense vector of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool TF >     // Transpose flag of the dense vector
inline const typename EnableIf< And< IsNumeric<ST1>, IsInvertible< typename DivTrait<ST1,ST2>::Type > >
                              , typename MultExprTrait< ST1, DVecScalarDivExpr<VT,ST2,TF> >::Type >::Type
   operator*( ST1 scalar, const DVecScalarDivExpr<VT,ST2,TF>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return vec.leftOperand() * ( scalar / vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a dense vector-scalar division expression
//        and a scalar value (\f$ \vec{a}=(\vec{b}/s1)/s2 \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector-scalar division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the division of a dense
// vector-scalar division expression and a scalar value.
*/
template< typename VT     // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the dense vector
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsNumeric<ST2>
                              , typename DivExprTrait<VT,typename MultTrait<ST1,ST2>::Type>::Type >::Type
   operator/( const DVecScalarDivExpr<VT,ST1,TF>& vec, ST2 scalar )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( scalar != ST2(0), "Division by zero detected" );

   typedef typename MultTrait<ST1,ST2>::Type         MultType;
   typedef typename DivExprTrait<VT,MultType>::Type  ReturnType;
   typedef typename ReturnType::RightOperand         ScalarType;

   if( IsMultExpr<ReturnType>::value ) {
      return ReturnType( vec.leftOperand(), ScalarType(1)/( vec.rightOperand() * scalar ) );
   }
   else {
      return ReturnType( vec.leftOperand(), vec.rightOperand() * scalar );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST, bool TF >
struct Size< DVecScalarDivExpr<VT,ST,TF> > : public Size<VT>
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
template< typename VT, typename ST, bool TF >
struct IsAligned< DVecScalarDivExpr<VT,ST,TF> > : public IsTrue< IsAligned<VT>::value >
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
template< typename VT, typename ST, bool TF >
struct IsPadded< DVecScalarDivExpr<VT,ST,TF> > : public IsTrue< IsPadded<VT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DVECSCALARMULTEXPRTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename ST2 >
struct DVecScalarMultExprTrait< DVecScalarDivExpr<VT,ST1,false>, ST2 >
{
 private:
   //**********************************************************************************************
   typedef typename DivTrait<ST2,ST1>::Type  ScalarType;
   //**********************************************************************************************

   //**********************************************************************************************
   enum { condition = IsInvertible<ScalarType>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename DVecScalarMultExprTrait<VT,ScalarType>::Type              T1;
   typedef DVecScalarMultExpr< DVecScalarDivExpr<VT,ST1,false>, ST2, false >  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsColumnVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SelectType<condition,T1,T2>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECSCALARMULTEXPRTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename ST2 >
struct TDVecScalarMultExprTrait< DVecScalarDivExpr<VT,ST1,true>, ST2 >
{
 private:
   //**********************************************************************************************
   typedef typename DivTrait<ST2,ST1>::Type  ScalarType;
   //**********************************************************************************************

   //**********************************************************************************************
   enum { condition = IsInvertible<ScalarType>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename DVecScalarMultExprTrait<VT,ScalarType>::Type            T1;
   typedef DVecScalarMultExpr< DVecScalarDivExpr<VT,ST1,true>, ST2, true >  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsRowVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SelectType<condition,T1,T2>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
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
template< typename VT, typename ST, bool TF, bool AF >
struct SubvectorExprTrait< DVecScalarDivExpr<VT,ST,TF>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename DivExprTrait< typename SubvectorExprTrait<const VT,AF>::Type, ST >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
