//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecScalarMultExpr.h
//  \brief Header file for the dense vector/scalar multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSCALARMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSCALARMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Expression.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/Assert.h>
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
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSCALARMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector-scalar multiplications.
// \ingroup dense_vector_expression
//
// The DVecScalarMultExpr class represents the compile time expression for multiplications between
// a dense vector and a scalar value.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar value
        , bool TF >    // Transpose flag
class DVecScalarMultExpr : public DenseVector< DVecScalarMultExpr<VT,ST,TF>, TF >
                         , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     RT;  //!< Result type of the dense vector expression.
   typedef typename VT::ReturnType     RN;  //!< Return type of the dense vector expression.
   typedef typename VT::ElementType    ET;  //!< Element type of the dense vector expression.
   typedef typename VT::CompositeType  CT;  //!< Composite type of the dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case the dense vector
       operand requires an intermediate evaluation, \a useAssign will be set to \a true and
       the multiplication expression will be evaluated via the \a assign function family.
       Otherwise \a useAssign will be set to \a false and the expression will be evaluated
       via the subscript operator. */
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
   typedef DVecScalarMultExpr<VT,ST,TF>                This;           //!< Type of this DVecScalarMultExpr instance.
   typedef typename MathTrait<RT,ST>::MultType         ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const typename MultExprTrait<RN,ST>::Type   ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DVecScalarMultExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  LeftOperand;

   //! Composite type of the right-hand side scalar value.
   typedef typename MathTrait< typename BaseElementType<VT>::Type, ST >::MultType  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT::vectorizable &&
                         IsSame<ET,RightOperand>::value &&
                         IntrinsicTrait<ET>::multiplication };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<VT>::value };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecScalarMultExpr class.
   //
   // \param vector The left-hand side dense vector of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   explicit inline DVecScalarMultExpr( const VT& vector, ST scalar )
      : vector_( vector )  // Left-hand side dense vector of the multiplication expression
      , scalar_( scalar )  // Right-hand side scalar of the multiplication expression
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
      return vector_[index] * scalar_;
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
      BLAZE_INTERNAL_ASSERT( index < vector_.size() , "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL, "Invalid vector access index" );
      const IntrinsicType xmm1( vector_.get( index ) );
      const IntrinsicType xmm2( set( scalar_ ) );
      return xmm1 * xmm2;
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
   LeftOperand  vector_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-scalar
   // multiplication expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the vector
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( ~lhs, rhs.vector_ );

      const size_t size( rhs.size() );
      for( size_t i=0UL; i<size; ++i )
         (~lhs)[i] *= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-scalar multiplication to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-scalar
   // multiplication expression to a sparse vector. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the vector
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( SparseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( ~lhs, rhs.vector_ );

      typename VT2::Iterator begin( (~lhs).begin() );
      const typename VT2::Iterator end( (~lhs).end() );

      for( ; begin!=end; ++begin )
         begin->value() *= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // scalar multiplication expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case the
   // vector operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      addAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
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
   /*!\brief Subtraction assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // scalar multiplication expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the vector
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      subAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
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
   /*!\brief Multiplication assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-scalar multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this operator can only be selected by the compiler in case the
   // vector operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      multAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL UNARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Unary minus operator for the negation of a dense vector (\f$ \vec{a} = -\vec{b} \f$).
// \ingroup dense_vector
//
// \param dv The dense vector to be negated.
// \return The negation of the vector.
//
// This operator represents the negation of a dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = -a;
   \endcode

// The operator returns an expression representing the negation of the given dense vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline const DVecScalarMultExpr<VT,typename BaseElementType<VT>::Type,TF>
   operator-( const DenseVector<VT,TF>& dv )
{
   typedef typename BaseElementType<VT>::Type  ElementType;
   return DVecScalarMultExpr<VT,ElementType,TF>( ~dv, ElementType(-1) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense vector and a scalar value
//        (\f$ \vec{a}=\vec{b}*s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator represents the multiplication between a dense vector and a scalar value:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = a * 1.25;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a T1::ElementType and \a T2. Both data types \a T1::ElementType and
// \a T2 have to be supported by the MathTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline const typename EnableIf< IsNumeric<T2>, typename MultExprTrait<T1,T2>::Type >::Type
   operator*( const DenseVector<T1,TF>& vec, T2 scalar )
{
   typedef typename MultExprTrait<T1,T2>::Type  Type;
   return Type( ~vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense vector
//        (\f$ \vec{a}=s*\vec{b} \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side vector for the multiplication.
// \return The scaled result vector.
//
// This operator represents the multiplication between a a scalar value and dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = 1.25 * a;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the involved data types \a T1::ElementType and \a T2. Both data types \a T1 and
// \a T2::ElementType have to be supported by the MathTrait class template. Note that this
// operator only works for scalar values of built-in data type.
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag
inline const typename EnableIf< IsNumeric<T1>, typename MultExprTrait<T1,T2>::Type >::Type
   operator*( T1 scalar, const DenseVector<T2,TF>& vec )
{
   typedef typename MultExprTrait<T1,T2>::Type  Type;
   return Type( ~vec, scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING UNARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unary minus operator for the negation of a dense vector-scalar multiplication
//        (\f$ \vec{a} = -(\vec{b} * s) \f$).
// \ingroup dense_vector
//
// \param dv The dense vector-scalar multiplication to be negated.
// \return The negation of the dense vector-scalar multiplication.
//
// This operator implements a performance optimized treatment of the negation of a dense vector-
// scalar multiplication expression.
*/
template< typename VT  // Type of the dense vector
        , typename ST  // Type of the scalar
        , bool TF >    // Transpose flag
inline const DVecScalarMultExpr<VT,ST,TF>
   operator-( const DVecScalarMultExpr<VT,ST,TF>& dv )
{
   return DVecScalarMultExpr<VT,ST,TF>( dv.leftOperand(), -dv.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar multiplication
//        expression and a scalar value (\f$ \vec{a}=(\vec{b}*s1)*s2 \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector-scalar multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// dense vector-scalar multiplication expression and a scalar value.
*/
template< typename VT     // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the dense vector
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsNumeric<ST2>
                              , typename MultExprTrait< DVecScalarMultExpr<VT,ST1,TF>, ST2 >::Type >::Type
   operator*( const DVecScalarMultExpr<VT,ST1,TF>& vec, ST2 scalar )
{
   return vec.leftOperand() * ( vec.rightOperand() * scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar multiplication
//        expression and a scalar value (\f$ \vec{a}=s2*(\vec{b}*s1) \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a dense vector-scalar multiplication expression.
*/
template< typename ST1  // Type of the left-hand side scalar
        , typename VT   // Type of the dense vector of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool TF >     // Transpose flag of the dense vector
inline const typename EnableIf< IsNumeric<ST1>
                              , typename MultExprTrait< ST1, DVecScalarMultExpr<VT,ST2,TF> >::Type >::Type
   operator*( ST1 scalar, const DVecScalarMultExpr<VT,ST2,TF>& vec )
{
   return vec.leftOperand() * ( scalar * vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a dense vector-scalar multiplication
//        expression by a scalar value (\f$ \vec{a}=(\vec{b}*s1)/s2 \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector-scalar multiplication.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the division of a
// dense vector-scalar multiplication expression by a scalar value.
*/
template< typename VT     // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the dense vector
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsFloatingPoint<typename MathTrait<ST1,ST2>::DivType>
                              , typename DivExprTrait< DVecScalarMultExpr<VT,ST1,TF>, ST2 >::Type >::Type
   operator/( const DVecScalarMultExpr<VT,ST1,TF>& vec, ST2 scalar )
{
   return vec.leftOperand() * ( vec.rightOperand() / scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar multiplication
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side dense vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense vector-scalar multiplication and a dense vector. It restructures the expression
// \f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST     // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the dense vectors
        , typename VT2 >  // Type of the right-hand side dense vector
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST,TF>, VT2 >::Type
   operator*( const DVecScalarMultExpr<VT1,ST,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   return ( lhs.leftOperand() * (~rhs) ) * lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector and a dense vector-
//        scalar multiplication expression (\f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense vector and a dense vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1   // Type of the left-hand side dense vector
        , bool TF        // Transpose flag of the dense vectors
        , typename VT2   // Type of the dense vector of the right-hand side expression
        , typename ST >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< VT1, DVecScalarMultExpr<VT2,ST,TF> >::Type
   operator*( const DenseVector<VT1,TF>& lhs, const DVecScalarMultExpr<VT2,ST,TF>& rhs )
{
   return ( (~lhs) * rhs.leftOperand() ) * rhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of two dense vector-scalar
//        multiplication expressions (\f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of
// two dense vector-scalar multiplication expressions. It restructures the expression
// \f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*(s1*s2) \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the dense vectors
        , typename VT2    // Type of the dense vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST1,TF>, DVecScalarMultExpr<VT2,ST2,TF> >::Type
   operator*( const DVecScalarMultExpr<VT1,ST1,TF>& lhs, const DVecScalarMultExpr<VT2,ST2,TF>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of a dense vector-scalar multiplication
//        expression and a dense vector (\f$ A=(\vec{b}*s1)*\vec{c}^T \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side dense vector.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the outer product of a
// dense vector-scalar multiplication and a dense vector. It restructures the expression
// \f$ A=(\vec{b}*s1)*\vec{c}^T \f$ to the expression \f$ A=(\vec{b}*\vec{c}^T)*s1 \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST     // Type of the scalar of the left-hand side expression
        , typename VT2 >  // Type of the right-hand side dense vector
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST,false>, VT2 >::Type
   operator*( const DVecScalarMultExpr<VT1,ST,false>& lhs, const DenseVector<VT2,true>& rhs )
{
   return ( lhs.leftOperand() * (~rhs) ) * lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of a dense vector and a dense vector-
//        scalar multiplication expression (\f$ A=\vec{b}*(\vec{c}^T*s1) \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the outer product of a
// dense vector and a dense vector-scalar multiplication. It restructures the expression
// \f$ A=\vec{b}*(\vec{c}^T*s1) \f$ to the expression \f$ A=(\vec{b}*\vec{c}^T)*s1 \f$.
*/
template< typename VT1   // Type of the left-hand side dense vector
        , typename VT2   // Type of the dense vector of the right-hand side expression
        , typename ST >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< VT1, DVecScalarMultExpr<VT2,ST,true> >::Type
   operator*( const DenseVector<VT1,false>& lhs, const DVecScalarMultExpr<VT2,ST,true>& rhs )
{
   return ( (~lhs) * rhs.leftOperand() ) * rhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of two a dense vector-scalar
//        multiplication expressions (\f$ A=(\vec{b}*s1)*(\vec{c}^T*s2) \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the outer product
// of two dense vector-scalar multiplications. It restructures the expression
// \f$ A=(\vec{b}*s1)*(\vec{c}^T*s2) \f$ to the expression \f$ A=(\vec{b}*\vec{c}^T)*(s1*s2) \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , typename VT2    // Type of the dense vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST1,false>, DVecScalarMultExpr<VT2,ST2,true> >::Type
   operator*( const DVecScalarMultExpr<VT1,ST1,false>& lhs, const DVecScalarMultExpr<VT2,ST2,true>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar multiplication
//        expression and a sparse vector (\f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side sparse vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense vector-scalar multiplication and a sparse vector. It restructures the expression
// \f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST     // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the vectors
        , typename VT2 >  // Type of the right-hand side sparse vector
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST,TF>, VT2 >::Type
   operator*( const DVecScalarMultExpr<VT1,ST,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   return ( lhs.leftOperand() * (~rhs) ) * lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse vector and a dense vector-
//        scalar multiplication expression (\f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// sparse vector and a dense vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1   // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag of the vectors
        , typename VT2   // Type of the dense vector of the right-hand side expression
        , typename ST >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< VT1, DVecScalarMultExpr<VT2,ST,TF> >::Type
   operator*( const SparseVector<VT1,TF>& lhs, const DVecScalarMultExpr<VT2,ST,TF>& rhs )
{
   return ( (~lhs) * rhs.leftOperand() ) * rhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar multiplication
//        expression and a sparse vector-scalar multiplication (\f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side sparse vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense vector-scalar multiplication and a sparse vector-scalar multiplication.
// It restructures the expression \f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$ to the
// expression \f$ \vec{a}=(\vec{b}*\vec{c})*(s1*s2) \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the vectors
        , typename VT2    // Type of the sparse vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar o the right-hand side expression
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST1,TF>, SVecScalarMultExpr<VT2,ST2,TF> >::Type
   operator*( const DVecScalarMultExpr<VT1,ST1,TF>& lhs, const SVecScalarMultExpr<VT2,ST2,TF>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse vector-scalar multiplication
//        expression and a dense vector-scalar multiplication (\f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector-scalar multiplication.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse vector-scalar multiplication and a dense vector-scalar multiplication.
// It restructures the expression \f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$ to the
// expression \f$ \vec{a}=(\vec{b}*\vec{c})*(s1*s2) \f$.
*/
template< typename VT1    // Type of the sparse vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool TF         // Transpose flag of the vectors
        , typename VT2    // Type of the dense vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar o the right-hand side expression
inline const typename MultExprTrait< SVecScalarMultExpr<VT1,ST1,TF>, DVecScalarMultExpr<VT2,ST2,TF> >::Type
   operator*( const SVecScalarMultExpr<VT1,ST1,TF>& lhs, const DVecScalarMultExpr<VT2,ST2,TF>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of a dense vector-scalar multiplication
//        expression and a sparse vector (\f$ A=(\vec{b}*s1)*\vec{c}^T \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side sparse vector.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the outer product of a
// dense vector-scalar multiplication and a sparse vector. It restructures the expression
// \f$ A=(\vec{b}*s1)*\vec{c}^T \f$ to the expression \f$ A=(\vec{b}*\vec{c}^T)*s1 \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST     // Type of the scalar of the left-hand side expression
        , typename VT2 >  // Type of the right-hand side sparse vector
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST,false>, VT2 >::Type
   operator*( const DVecScalarMultExpr<VT1,ST,false>& lhs, const SparseVector<VT2,true>& rhs )
{
   return ( lhs.leftOperand() * (~rhs) ) * lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of a sparse vector and a dense vector-
//        scalar multiplication expression (\f$ A=\vec{b}*(\vec{c}^T*s1) \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse vector.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the outer product of a
// sparse vector and a dense vector-scalar multiplication. It restructures the expression
// \f$ A=\vec{b}*(\vec{c}^T*s1) \f$ to the expression \f$ A=(\vec{b}*\vec{c}^T)*s1 \f$.
*/
template< typename VT1   // Type of the left-hand side sparse vector
        , typename VT2   // Type of the dense vector of the right-hand side expression
        , typename ST >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< VT1, DVecScalarMultExpr<VT2,ST,true> >::Type
   operator*( const SparseVector<VT1,false>& lhs, const DVecScalarMultExpr<VT2,ST,true>& rhs )
{
   return ( (~lhs) * rhs.leftOperand() ) * rhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of a dense vector-scalar multiplication
//        expression and a sparse vector-scalar multiplication (\f$ A=(\vec{b}*s1)*(\vec{c}^T*s2) \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side dense vector-scalar multiplication.
// \param rhs The right-hand side sparse vector-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the outer product
// of a dense vector-scalar multiplication and a sparse vector-scalar multiplication.
// It restructures the expression \f$ A=(\vec{b}*s1)*(\vec{c}^T*s2) \f$ to the
// expression \f$ A=(\vec{b}*\vec{c}^T)*(s1*s2) \f$.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , typename VT2    // Type of the sparse vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar o the right-hand side expression
inline const typename MultExprTrait< DVecScalarMultExpr<VT1,ST1,false>, SVecScalarMultExpr<VT2,ST2,true> >::Type
   operator*( const DVecScalarMultExpr<VT1,ST1,false>& lhs, const SVecScalarMultExpr<VT2,ST2,true>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the outer product of a sparse vector-scalar multiplication
//        expression and a dense vector-scalar multiplication (\f$ A=(\vec{b}*s1)*(\vec{c}^T*s2) \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse vector-scalar multiplication.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse vector-scalar multiplication and a dense vector-scalar multiplication.
// It restructures the expression \f$ A=(\vec{b}*s1)*(\vec{c}^T*s2) \f$ to the
// expression \f$ A=(\vec{b}*\vec{c}^T)*(s1*s2) \f$.
*/
template< typename VT1    // Type of the sparse vector of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , typename VT2    // Type of the dense vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar o the right-hand side expression
inline const typename MultExprTrait< SVecScalarMultExpr<VT1,ST1,false>, DVecScalarMultExpr<VT2,ST2,true> >::Type
   operator*( const SVecScalarMultExpr<VT1,ST1,false>& lhs, const DVecScalarMultExpr<VT2,ST2,true>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix and a dense
//        vector-scalar multiplication expression (\f$ \vec{a}=B*(\vec{c}*s1) \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense matrix.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense matrix and a dense vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=B*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order of the left-hand side dense matrix
        , typename VT    // Type of the dense vector of the right-hand side expression
        , typename ST >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< MT, DVecScalarMultExpr<VT,ST,false> >::Type
   operator*( const DenseMatrix<MT,SO>& mat, const DVecScalarMultExpr<VT,ST,false>& vec )
{
   return ( (~mat) * vec.leftOperand() ) * vec.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose dense vector-scalar
//        multiplication expression and a dense matrix (\f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side transpose dense vector-scalar multiplication.
// \param rhs The right-hand side dense matrix.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// transpose dense vector-scalar multiplication and a dense matrix. It restructures the
// expression \f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$ to the expression \f$ \vec{a}^T=(\vec{b}^T*C)*s1 \f$.
*/
template< typename VT  // Type of the dense vector of the left-hand side expression
        , typename ST  // Type of the scalar of the left-hand side expression
        , typename MT  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order of the right-hand side dense matrix
inline const typename MultExprTrait< DVecScalarMultExpr<VT,ST,true>, MT >::Type
   operator*( const DVecScalarMultExpr<VT,ST,true>& vec, const DenseMatrix<MT,SO>& mat )
{
   return ( vec.leftOperand() * (~mat) ) * vec.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse matrix and a dense
//        vector-scalar multiplication expression (\f$ \vec{a}=B*(\vec{c}*s1) \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse matrix.
// \param rhs The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// sparse matrix and a dense vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=B*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order of the left-hand side sparse matrix
        , typename VT    // Type of the dense vector of the right-hand side expression
        , typename ST >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< MT, DVecScalarMultExpr<VT,ST,false> >::Type
   operator*( const SparseMatrix<MT,SO>& mat, const DVecScalarMultExpr<VT,ST,false>& vec )
{
   return ( (~mat) * vec.leftOperand() ) * vec.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose dense vector-scalar
//        multiplication expression and a sparse matrix (\f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side transpose dense vector-scalar multiplication.
// \param rhs The right-hand side sparse matrix.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// transpose dense vector-scalar multiplication and a sparse matrix. It restructures the
// expression \f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$ to the expression \f$ \vec{a}^T=(\vec{b}^T*C)*s1 \f$.
*/
template< typename VT  // Type of the dense vector of the left-hand side expression
        , typename ST  // Type of the scalar of the left-hand side expression
        , typename MT  // Type of the right-hand side sparse matrix
        , bool SO >    // Storage order of the right-hand side sparse matrix
inline const typename MultExprTrait< DVecScalarMultExpr<VT,ST,true>, MT >::Type
   operator*( const DVecScalarMultExpr<VT,ST,true>& vec, const SparseMatrix<MT,SO>& mat )
{
   return ( vec.leftOperand() * (~mat) ) * vec.rightOperand();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DVECSCALARMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename ST2 >
struct DVecScalarMultTrait< DVecScalarMultExpr<VT,ST1,false>, ST2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DVecScalarMultTrait<VT,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECSCALARMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename ST2 >
struct TDVecScalarMultTrait< DVecScalarMultExpr<VT,ST1,true>, ST2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDVecScalarMultTrait<VT,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DVECSCALARDIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename ST2 >
struct DVecScalarDivTrait< DVecScalarMultExpr<VT,ST1,false>, ST2 >
{
 private:
   //**********************************************************************************************
   enum { condition = IsFloatingPoint<typename MathTrait<ST1,ST2>::DivType>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename DVecScalarMultTrait<VT,typename MathTrait<ST1,ST2>::DivType>::Type  T1;
   typedef typename DVecScalarDivTrait<VT,typename MathTrait<ST1,ST2>::DivType>::Type   T2;
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
//  TDVECSCALARDIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename ST2 >
struct TDVecScalarDivTrait< DVecScalarMultExpr<VT,ST1,true>, ST2 >
{
 private:
   //**********************************************************************************************
   enum { condition = IsFloatingPoint<typename MathTrait<ST1,ST2>::DivType>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename TDVecScalarMultTrait<VT,typename MathTrait<ST1,ST2>::DivType>::Type  T1;
   typedef typename TDVecScalarDivTrait<VT,typename MathTrait<ST1,ST2>::DivType>::Type   T2;
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




//=================================================================================================
//
//  DVECDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST, typename VT2 >
struct DVecDVecMultTrait< DVecScalarMultExpr<VT1,ST,false>, VT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename DVecDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct DVecDVecMultTrait< VT1, DVecScalarMultExpr<VT2,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename DVecDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct DVecDVecMultTrait< DVecScalarMultExpr<VT1,ST1,false>, DVecScalarMultExpr<VT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DVecScalarMultTrait<typename DVecDVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DVECTDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST, typename VT2 >
struct DVecTDVecMultTrait< DVecScalarMultExpr<VT1,ST,false>, VT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DVecTDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct DVecTDVecMultTrait< VT1, DVecScalarMultExpr<VT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DVecTDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct DVecTDVecMultTrait< DVecScalarMultExpr<VT1,ST1,false>, DVecScalarMultExpr<VT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename DVecTDVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECTDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST, typename VT2 >
struct TDVecTDVecMultTrait< DVecScalarMultExpr<VT1,ST,true>, VT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecTDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct TDVecTDVecMultTrait< VT1, DVecScalarMultExpr<VT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecTDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct TDVecTDVecMultTrait< DVecScalarMultExpr<VT1,ST1,true>, DVecScalarMultExpr<VT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value && IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDVecScalarMultTrait<typename TDVecTDVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DVECSVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct DVecSVecMultTrait< DVecScalarMultExpr<VT1,ST,false>, VT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value  && !IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename SVecScalarMultTrait<typename DVecSVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct DVecSVecMultTrait< DVecScalarMultExpr<VT1,ST1,false>, SVecScalarMultExpr<VT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value  && !IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SVecScalarMultTrait<typename DVecSVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DVECTSVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST, typename VT2 >
struct DVecTSVecMultTrait< DVecScalarMultExpr<VT1,ST,false>, VT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value  && !IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST>::value
                              , typename SMatScalarMultTrait<typename DVecTSVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct DVecTSVecMultTrait< DVecScalarMultExpr<VT1,ST1,false>, SVecScalarMultExpr<VT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value  && !IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SMatScalarMultTrait<typename DVecTSVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECTSVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST, typename VT2 >
struct TDVecTSVecMultTrait< DVecScalarMultExpr<VT1,ST,true>, VT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value  && IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename TSVecScalarMultTrait<typename TDVecTSVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct TDVecTSVecMultTrait< DVecScalarMultExpr<VT1,ST1,true>, SVecScalarMultExpr<VT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT1>::value  && IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TSVecScalarMultTrait<typename TDVecTSVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SVECDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct SVecDVecMultTrait< VT1, DVecScalarMultExpr<VT2,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename SVecScalarMultTrait<typename SVecDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct SVecDVecMultTrait< SVecScalarMultExpr<VT1,ST1,false>, DVecScalarMultExpr<VT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && !IsTransposeVector<VT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SVecScalarMultTrait<typename SVecDVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SVECTDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct SVecTDVecMultTrait< VT1, DVecScalarMultExpr<VT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST>::value
                              , typename TSMatScalarMultTrait<typename SVecTDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct SVecTDVecMultTrait< SVecScalarMultExpr<VT1,ST1,false>, DVecScalarMultExpr<VT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT1>::value && !IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsTransposeVector<VT2>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TSMatScalarMultTrait<typename SVecTDVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TSVECTDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename ST >
struct TSVecTDVecMultTrait< VT1, DVecScalarMultExpr<VT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT1>::value && IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST>::value
                              , typename TSVecScalarMultTrait<typename TSVecTDVecMultTrait<VT1,VT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename ST1, typename VT2, typename ST2 >
struct TSVecTDVecMultTrait< SVecScalarMultExpr<VT1,ST1,true>, DVecScalarMultExpr<VT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT1>::value && IsTransposeVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsTransposeVector<VT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TSVecScalarMultTrait<typename TSVecTDVecMultTrait<VT1,VT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT, typename ST >
struct DMatDVecMultTrait< MT, DVecScalarMultExpr<VT,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value   &&
                                IsDenseVector<VT>::value && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename DMatDVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT, typename ST >
struct TDMatDVecMultTrait< MT, DVecScalarMultExpr<VT,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsDenseVector<VT>::value && !IsTransposeVector<VT>::value  &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename TDMatDVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TDVecDMatMultTrait< DVecScalarMultExpr<VT,ST,true>, MT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value  &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecDMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TDVecTDMatMultTrait< DVecScalarMultExpr<VT,ST,true>, MT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value   &&
                                IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecTDMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SMATDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT, typename ST >
struct SMatDVecMultTrait< MT, DVecScalarMultExpr<VT,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT>::value && IsRowMajorMatrix<MT>::value   &&
                                IsDenseVector<VT>::value  && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename SMatDVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TSMATDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT, typename ST >
struct TSMatDVecMultTrait< MT, DVecScalarMultExpr<VT,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsDenseVector<VT>::value  && !IsTransposeVector<VT>::value  &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename TSMatDVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECSMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TDVecSMatMultTrait< DVecScalarMultExpr<VT,ST,true>, MT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value  && IsTransposeVector<VT>::value &&
                                IsSparseMatrix<MT>::value && IsRowMajorMatrix<MT>::value  &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecSMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECTSMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TDVecTSMatMultTrait< DVecScalarMultExpr<VT,ST,true>, MT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value  && IsTransposeVector<VT>::value   &&
                                IsSparseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecTSMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
