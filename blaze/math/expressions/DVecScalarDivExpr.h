//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecScalarDivExpr.h
//  \brief Header file for the dense vector/scalar division expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSCALARDIVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSCALARDIVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Expression.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSCALARDIVEXPRHELPER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Helper class for divisions of a dense vector by a scalar.
// \ingroup dense_vector_expression
//
// The DVecScalarDivExprHelper class is an auxiliary class to define the return type of the
// division between a dense vector and a scalar value.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar value
        , bool TF >    // Transpose flag
struct DVecScalarDivExprHelper
{
 public:
   //**Type definitions****************************************************************************
   //! Scalar type for the instantiation of the resulting expression object.
   typedef typename DivTrait< typename BaseElementType<VT>::Type, ST >::Type  ScalarType;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation of the dense vector/scalar division return type.
   enum { value = IsFloatingPoint<ScalarType>::value };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Resulting type of the division between the given dense vector and scalar value.
   typedef typename SelectType< value,
                                DVecScalarMultExpr<VT,ScalarType,TF>,
                                DVecScalarDivExpr<VT,ScalarType,TF> >::Type  Type;
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ScalarType );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




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
                        , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     RT;  //!< Result type of the dense vector expression.
   typedef typename VT::ReturnType     RN;  //!< Return type of the dense vector expression.
   typedef typename VT::CompositeType  CT;  //!< Composite type of the dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case either the dense vector
       operand requires an intermediate evaluation, \a useAssign will be set to \a true and the
       multiplication expression will be evaluated via the \a assign function family. Otherwise
       \a useAssign will be set to \a false and the expression will be evaluated via the subscript
       operator. */
   enum { useAssign = !IsReference<CT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct UseAssign {
      enum { value = useAssign };
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
   typedef const typename DivExprTrait<RN,ST>::Type    ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DVecScalarDivExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  LeftOperand;

   //! Composite type of the right-hand side scalar value.
   typedef typename DivTrait< typename BaseElementType<VT>::Type, ST >::Type  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<VT>::value };
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
   // principle, this operator can only be selected by the compiler in case the vector
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( ~lhs, rhs.vector_ );
      (~lhs) /= rhs.scalar_;
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
   // principle, this operator can only be selected by the compiler in case the vector
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( SparseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
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
   // SFINAE principle, this operator can only be selected by the compiler in case the vector
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      addAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
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
   // principle, this operator can only be selected by the compiler in case the vector operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      subAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
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
   // of the SFINAE principle, this operator can only be selected by the compiler in case the
   // vector operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      multAssign( DenseVector<VT2,TF>& lhs, const DVecScalarDivExpr& rhs )
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

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
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
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline const typename EnableIf< IsNumeric<T2>,
                                typename DVecScalarDivExprHelper<T1,T2,TF>::Type >::Type
   operator/( const DenseVector<T1,TF>& vec, T2 scalar )
{
   BLAZE_USER_ASSERT( scalar != T2(0), "Division by zero detected" );

   typedef DVecScalarDivExprHelper<T1,T2,TF>  Helper;
   typedef typename Helper::ScalarType        ScalarType;

   if( Helper::value ) {
      return typename Helper::Type( ~vec, ScalarType(1)/ScalarType(scalar) );
   }
   else {
      return typename Helper::Type( ~vec, scalar );
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
inline const typename EnableIf< IsFloatingPoint< typename DivTrait<ST2,ST1>::Type >
                              , typename MultExprTrait< DVecScalarDivExpr<VT,ST1,TF>, ST2 >::Type >::Type
   operator*( const DVecScalarDivExpr<VT,ST1,TF>& vec, ST2 scalar )
{
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
inline const typename EnableIf< IsFloatingPoint< typename DivTrait<ST1,ST2>::Type >
                              , typename MultExprTrait< ST1, DVecScalarDivExpr<VT,ST2,TF> >::Type >::Type
   operator*( ST1 scalar, const DVecScalarDivExpr<VT,ST2,TF>& vec )
{
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
                              , typename DVecScalarDivExprHelper<VT,typename MultTrait<ST1,ST2>::Type,TF>::Type >::Type
   operator/( const DVecScalarDivExpr<VT,ST1,TF>& vec, ST2 scalar )
{
   BLAZE_USER_ASSERT( scalar != ST2(0), "Division by zero detected" );

   typedef typename MultTrait<ST1,ST2>::Type        MultType;
   typedef DVecScalarDivExprHelper<VT,MultType,TF>  Helper;

   if( Helper::value ) {
      return typename Helper::Type( vec.leftOperand(), MultType(1)/( vec.rightOperand() * scalar ) );
   }
   else {
      return typename Helper::Type( vec.leftOperand(), vec.rightOperand() * scalar );
   }
}
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
   enum { condition = IsFloatingPoint<typename DivTrait<ST1,ST2>::Type>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename DVecScalarMultExprTrait<VT,typename DivTrait<ST1,ST2>::Type>::Type  T1;
   typedef DVecScalarMultExpr< DVecScalarDivExpr<VT,ST1,false>, ST2, false >            T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && !IsTransposeVector<VT>::value &&
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
   enum { condition = IsFloatingPoint<typename DivTrait<ST1,ST2>::Type>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename DVecScalarMultExprTrait<VT,typename DivTrait<ST1,ST2>::Type>::Type  T1;
   typedef DVecScalarMultExpr< DVecScalarDivExpr<VT,ST1,true>, ST2, true >              T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SelectType<condition,T1,T2>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
