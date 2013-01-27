//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecSVecAddExpr.h
//  \brief Header file for the dense vector/sparse vector addition expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSVECADDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSVECADDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSVECADDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector-sparse vector additions.
// \ingroup dense_vector_expression
//
// The DVecSVecAddExpr class represents the compile time expression for additions between
// a dense vector and a sparse vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
class DVecSVecAddExpr : public DenseVector< DVecSVecAddExpr<VT1,VT2,TF>, TF >
                      , private Expression
                      , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT1::ResultType     RT1;  //!< Result type of the left-hand side dense vector expression.
   typedef typename VT2::ResultType     RT2;  //!< Result type of the right-hand side sparse vector expression.
   typedef typename VT1::ReturnType     RN1;  //!< Return type of the left-hand side dense vector expression.
   typedef typename VT2::ReturnType     RN2;  //!< Return type of the right-hand side sparse vector expression.
   typedef typename VT1::CompositeType  CT1;  //!< Composite type of the left-hand side dense vector expression.
   typedef typename VT2::CompositeType  CT2;  //!< Composite type of the right-hand side sparse vector expression.
   typedef typename VT1::TransposeType  TT1;  //!< Transpose type of the left-hand side dense vector expression.
   typedef typename VT2::TransposeType  TT2;  //!< Transpose type of the right-hand side sparse vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   enum { returnExpr = !IsTemporary<RN1>::value && !IsTemporary<RN2>::value };

   //! Expression return type for the subscript operator.
   typedef typename AddExprTrait<RN1,RN2>::Type  ExprReturnType;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DVecSVecAddExpr<VT1,VT2,TF>         This;           //!< Type of this DVecSVecAddExpr instance.
   typedef typename AddTrait<RT1,RT2>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef const ResultType  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT1>::value, const VT1, const VT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse vector expression.
   typedef typename SelectType< IsExpression<VT2>::value, const VT2, const VT2& >::Type  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = ( IsComputation<VT1>::value && CanAlias<VT1>::value ) };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecSVecAddExpr class.
   //
   // \param lhs The left-hand side dense vector operand of the addition expression.
   // \param rhs The right-hand side sparse vector operand of the addition expression.
   */
   explicit inline DVecSVecAddExpr( const VT1& lhs, const VT2& rhs )
      : lhs_( lhs )  // Left-hand side dense vector of the addition expression
      , rhs_( rhs )  // Right-hand side sparse vector of the addition expression
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
      return lhs_[index] + rhs_[index];
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
   /*!\brief Returns the right-hand side sparse vector operand.
   //
   // \return The right-hand side sparse vector operand.
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
      return ( IsComputation<VT1>::value && CanAlias<VT1>::value && lhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense vector of the addition expression.
   RightOperand rhs_;  //!< Right-hand side sparse vector of the addition expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-sparse
   // vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign   ( ~lhs, rhs.lhs_ );
      addAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-sparse vector addition to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-sparse
   // vector addition expression to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

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
   /*!\brief Addition assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( ~lhs, rhs.lhs_ );
      addAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( ~lhs, rhs.lhs_ );
      subAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2 );
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
/*!\brief Addition operator for the addition of a dense vector and a sparse vector
//        (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the vector addition.
// \param rhs The right-hand side sparse vector for the vector addition.
// \return The sum of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the addition of a dense vector and a sparse vector:

   \code
   blaze::DynamicVector<double> a, c;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization
   c = a + b;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved vector element types \a T1::ElementType and \a T2::ElementType.
// Both vector types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side sparse vector
        , bool TF >    // Transpose flag
inline const DVecSVecAddExpr<T1,T2,TF>
   operator+( const DenseVector<T1,TF>& lhs, const SparseVector<T2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   return DVecSVecAddExpr<T1,T2,TF>( ~lhs, ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a sparse vector and a dense vector
//        (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector for the vector addition.
// \param rhs The right-hand side dense vector for the vector addition.
// \return The sum of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the addition of a sparse vector and a dense vector:

   \code
   blaze::CompressedVector<double> a;
   blaze::DynamicVector<double> b, c;
   // ... Resizing and initialization
   c = a + b;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved vector element types \a T1::ElementType and \a T2::ElementType.
// Both vector types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1  // Type of the left-hand side sparse vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag
inline const DVecSVecAddExpr<T2,T1,TF>
   operator+( const SparseVector<T1,TF>& lhs, const DenseVector<T2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   return DVecSVecAddExpr<T2,T1,TF>( ~rhs, ~lhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition operator for the addition of a dense vector-sparse vector addition
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}+\vec{c})+\vec{d} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector-sparse vector addition.
// \param rhs The right-hand side dense vector.
// \return The sum of the two vectors.
//
// This operator implements a performance optimized treatment of the addition of a dense
// vector-sparse vector addition expression to a dense vector.
*/
template< typename T1    // Type of the dense vector of the left-hand side expression
        , typename T2    // Type of the sparse vector of the left-hand side expression
        , bool TF        // Transpose flag of the left-hand side expression
        , typename T3 >  // Type of the right-hand side dense vector
inline const typename AddExprTrait< DVecSVecAddExpr<T1,T2,TF>, T3 >::Type
   operator+( const DVecSVecAddExpr<T1,T2,TF>& lhs, const DenseVector<T3,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() + (~rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction operator for the subtraction of a dense vector-sparse vector addition
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}+\vec{c})-\vec{d} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector-sparse vector addition.
// \param rhs The right-hand side dense vector.
// \return The difference of the two vectors.
//
// This operator implements a performance optimized treatment of the subtraction of a
// dense vector-sparse vector addition expression and a dense vector.
*/
template< typename T1    // Type of the dense vector of the left-hand side expression
        , typename T2    // Type of the sparse vector of the left-hand side expression
        , bool TF        // Transpose flag of the left-hand side expression
        , typename T3 >  // Type of the right-hand side dense vector
inline const typename SubExprTrait< DVecSVecAddExpr<T1,T2,TF>, T3 >::Type
   operator-( const DVecSVecAddExpr<T1,T2,TF>& lhs, const DenseVector<T3,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() - (~rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct DVecDVecAddExprTrait< DVecSVecAddExpr<VT1,VT2,false>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseVector<VT1>::value  && !IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && !IsTransposeVector<VT3>::value
                              , typename DVecSVecAddExprTrait< typename DVecDVecAddExprTrait<VT1,VT3>::Type, VT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct TDVecTDVecAddExprTrait< DVecSVecAddExpr<VT1,VT2,true>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseVector<VT1>::value  && IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && IsTransposeVector<VT3>::value
                              , typename TDVecTSVecAddExprTrait< typename TDVecTDVecAddExprTrait<VT1,VT3>::Type, VT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct DVecDVecSubExprTrait< DVecSVecAddExpr<VT1,VT2,false>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseVector<VT1>::value  && !IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && !IsTransposeVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && !IsTransposeVector<VT3>::value
                              , typename DVecSVecAddExprTrait< typename DVecDVecSubExprTrait<VT1,VT3>::Type, VT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct TDVecTDVecSubExprTrait< DVecSVecAddExpr<VT1,VT2,true>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseVector<VT1>::value  && IsTransposeVector<VT1>::value &&
                                IsSparseVector<VT2>::value && IsTransposeVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && IsTransposeVector<VT3>::value
                              , typename TDVecTSVecAddExprTrait< typename TDVecTDVecSubExprTrait<VT1,VT3>::Type, VT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
