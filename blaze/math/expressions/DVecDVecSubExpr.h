//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecDVecSubExpr.h
//  \brief Header file for the dense vector/dense vector subtraction expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECDVECSUBEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECDVECSUBEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Expression.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECDVECSUBEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector-dense vector subtractions.
// \ingroup dense_vector_expression
//
// The DVecDVecSubExpr class represents the compile time expression for subtractions between
// dense vectors.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
class DVecDVecSubExpr : public DenseVector< DVecDVecSubExpr<VT1,VT2,TF>, TF >
                      , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT1::ResultType     RT1;  //!< Result type of the left-hand side dense vector expression.
   typedef typename VT2::ResultType     RT2;  //!< Result type of the right-hand side dense vector expression.
   typedef typename VT1::ReturnType     RN1;  //!< Return type of the left-hand side dense vector expression.
   typedef typename VT2::ReturnType     RN2;  //!< Return type of the right-hand side dense vector expression.
   typedef typename VT1::CompositeType  CT1;  //!< Composite type of the left-hand side dense vector expression.
   typedef typename VT2::CompositeType  CT2;  //!< Composite type of the right-hand side dense vector expression.
   typedef typename VT1::ElementType    ET1;  //!< Element type of the left-hand side dense vector expression.
   typedef typename VT2::ElementType    ET2;  //!< Element type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the subtraction expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for the
       evaluation strategy of the subtraction expression. In case either of the two dense vector
       operands requires an intermediate evaluation, \a useAssign will be set to \a true and
       the subtraction expression will be evaluated via the \a assign function family. Otherwise
       \a useAssign will be set to \a false and the expression will be evaluated via the subscript
       operator. */
   enum { useAssign = ( !IsReference<CT1>::value || !IsReference<CT2>::value ) };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct UseAssign {
      enum { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DVecDVecSubExpr<VT1,VT2,TF>                 This;           //!< Type of this DVecDVecSubExpr instance.
   typedef typename MathTrait<RT1,RT2>::SubType        ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const typename SubExprTrait<RN1,RN2>::Type  ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DVecDVecSubExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT1>::value, const VT1, const VT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT2>::value, const VT2, const VT2& >::Type  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT1::vectorizable && VT2::vectorizable &&
                         IsSame<ET1,ET2>::value &&
                         IntrinsicTrait<ET1>::subtraction };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = ( IsExpression<VT1>::value && CanAlias<VT1>::value ) ||
                     ( IsExpression<VT2>::value && CanAlias<VT2>::value ) };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecDVecSubExpr class.
   //
   // \param lhs The left-hand side operand of the subtraction expression.
   // \param rhs The right-hand side operand of the subtraction expression.
   */
   explicit inline DVecDVecSubExpr( const VT1& lhs, const VT2& rhs )
      : lhs_( lhs )  // Left-hand side dense vector of the subtraction expression
      , rhs_( rhs )  // Right-hand side dense vector of the subtraction expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.size() == rhs.size(), "Invalid vector sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < lhs_.size(), "Invalid vector access index" );
      return lhs_[index] - rhs_[index];
   }
   //**********************************************************************************************

   //**Get function********************************************************************************
   /*!\brief Access to the intrinsic elements of the vector.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   inline IntrinsicType get( size_t index ) const {
      typedef IntrinsicTrait<ElementType>  IT;
      BLAZE_INTERNAL_ASSERT( index < lhs_.size()    , "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL, "Invalid vector access index" );
      const IntrinsicType xmm1( lhs_.get( index ) );
      const IntrinsicType xmm2( rhs_.get( index ) );
      return xmm1 - xmm2;
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return lhs_.size();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const {
      return rhs_;
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
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense vector of the subtraction expression.
   RightOperand rhs_;  //!< Right-hand side dense vector of the subtraction expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector subtraction expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT> >::Type
      assign( DenseVector<VT,TF>& lhs, const DVecDVecSubExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign   ( ~lhs, rhs.lhs_ );
      subAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector subtraction to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector subtraction expression to a sparse vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT> >::Type
      assign( SparseVector<VT,TF>& lhs, const DVecDVecSubExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // dense vector subtraction expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this operator can only be selected by the compiler in case either
   // of the operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT> >::Type
      addAssign( DenseVector<VT,TF>& lhs, const DVecDVecSubExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( ~lhs, rhs.lhs_ );
      subAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // dense vector subtraction expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case either of
   // the operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT> >::Type
      subAssign( DenseVector<VT,TF>& lhs, const DVecDVecSubExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( ~lhs, rhs.lhs_ );
      addAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-dense vector subtraction expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler
   // in case either of the operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT> >::Type
      multAssign( DenseVector<VT,TF>& lhs, const DVecDVecSubExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      multAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT1, TF );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );
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
/*!\brief Subtraction operator for the subtraction of two dense vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the vector subtraction.
// \param rhs The right-hand side dense vector to be subtracted from the vector.
// \return The difference of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the subtraction of two dense vectors:

   \code
   blaze::DynamicVector<double> a, b, c;
   // ... Resizing and initialization
   c = a - b;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved vector element types \a T1::ElementType and \a T2::ElementType.
// Both vector types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the MathTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag
inline const DVecDVecSubExpr<T1,T2,TF>
   operator-( const DenseVector<T1,TF>& lhs, const DenseVector<T2,TF>& rhs )
{
   if( (~lhs).size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   return DVecDVecSubExpr<T1,T2,TF>( ~lhs, ~rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
