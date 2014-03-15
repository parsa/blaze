//=================================================================================================
/*!
//  \file blaze/math/expressions/TSVecDMatMultExpr.h
//  \brief Header file for the transpose sparse vector/dense matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TSVECDMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSVECDMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TSVECDMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose sparse vector-dense matrix multiplications.
// \ingroup dense_vector_expression
//
// The TSVecDMatMultExpr class represents the compile time expression for multiplications
// between transpose sparse vectors and row-major dense matrices.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side dense matrix
class TSVecDMatMultExpr : public DenseVector< TSVecDMatMultExpr<VT,MT>, true >
                        , private TVecMatMultExpr
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     VRT;  //!< Result type of the left-hand side sparse vector expression.
   typedef typename MT::ResultType     MRT;  //!< Result type of the right-hand side dense matrix expression.
   typedef typename VRT::ElementType   VET;  //!< Element type of the left-hand side sparse vector expression.
   typedef typename MRT::ElementType   MET;  //!< Element type of the right-hand side dense matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the left-hand side sparse vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the right-hand side dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense vector expression.
   enum { evaluateVector = IsComputation<VT>::value || RequiresEvaluation<VT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   enum { evaluateMatrix = ( IsComputation<MT>::value && IsSame<MET,VET>::value &&
                             IsBlasCompatible<MET>::value ) || RequiresEvaluation<MT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the either the vector or the matrix operand require an intermediate evaluation,
       the nested \value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseSMPAssignKernel {
      enum { value = evaluateVector || evaluateMatrix };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the matrix type and the two involved vector types are suited for a vectorized
       computation of the vector/matrix multiplication, the nested \value will be set to 1,
       otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseVectorizedKernel {
      enum { value = !UseSMPAssignKernel<T1,T2,T3>::value &&
                     T1::vectorizable && T3::vectorizable &&
                     IsSame<typename T1::ElementType,typename T2::ElementType>::value &&
                     IsSame<typename T1::ElementType,typename T3::ElementType>::value &&
                     IntrinsicTrait<typename T1::ElementType>::addition &&
                     IntrinsicTrait<typename T1::ElementType>::multiplication };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case a vectorized computation of the vector/matrix multiplication is not possible, but
       a loop-unrolled computation is feasible, the nested \value will be set to 1, otherwise it
       will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseOptimizedKernel {
      enum { value = !UseSMPAssignKernel<T1,T2,T3>::value &&
                     !UseVectorizedKernel<T1,T2,T3>::value &&
                     !IsResizable<typename T1::ElementType>::value &&
                     !IsResizable<VET>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case neither a vectorized nor optimized computation is possible, the nested \value will
       be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseDefaultKernel {
      enum { value = !UseSMPAssignKernel<T1,T2,T3>::value &&
                     !UseVectorizedKernel<T1,T2,T3>::value &&
                     !UseOptimizedKernel<T1,T2,T3>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TSVecDMatMultExpr<VT,MT>                    This;           //!< Type of this TSVecDMatMultExpr instance.
   typedef typename MultTrait<VRT,MRT>::Type           ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const ElementType                           ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                            CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side sparse vector operand.
   typedef typename SelectType< evaluateVector, const VRT, VCT >::Type  LT;

   //! Type for the assignment of the right-hand side dense matrix operand.
   typedef typename SelectType< evaluateMatrix, const MRT, MCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = MT::vectorizable &&
                         IsSame<VET,MET>::value &&
                         IntrinsicTrait<VET>::addition &&
                         IntrinsicTrait<VET>::multiplication };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateVector && !evaluateMatrix };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TSVecDMatMultExpr class.
   //
   // \param vec The left-hand side sparse vector operand of the multiplication expression.
   // \param mat The right-hand side dense matrix operand of the multiplication expression.
   */
   explicit inline TSVecDMatMultExpr( const VT& vec, const MT& mat )
      : vec_( vec )  // Left-hand side sparse vector of the multiplication expression
      , mat_( mat )  // Right-hand side dense matrix of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( vec_.size() == mat_.rows(), "Invalid vector and matrix sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < mat_.columns(), "Invalid vector access index" );

      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      LT x( vec_ );  // Evaluation of the left-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == vec_.size(), "Invalid vector size" );

      const ConstIterator end( x.end() );
      ConstIterator element( x.begin() );
      ElementType res;

      if( element != end ) {
         res = element->value() * mat_( element->index(), index );
         ++element;
         for( ; element!=end; ++element )
            res += element->value() * mat_( element->index(), index );
      }
      else {
         reset( res );
      }

      return res;
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return mat_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse vector operand.
   //
   // \return The left-hand side sparse vector operand.
   */
   inline LeftOperand leftOperand() const {
      return vec_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense matrix operand.
   //
   // \return The right-hand side dense matrix operand.
   */
   inline RightOperand rightOperand() const {
      return mat_;
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
      return vec_.isAliased( alias ) || mat_.isAliased( alias );
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
      return vec_.isAliased( alias ) || mat_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const {
      return mat_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return ( size() > SMP_TSVECDMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vec_;  //!< Left-hand side sparse vector of the multiplication expression.
   RightOperand mat_;  //!< Right-hand side dense matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse vector-dense matrix multiplication to a dense vector
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse vector-
   // dense matrix multiplication expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT2,true>& lhs, const TSVecDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) {
         reset( ~lhs );
         return;
      }

      // Evaluation of the right-hand side dense matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse vector-dense matrix multiplication
      TSVecDMatMultExpr::selectAssignKernel( ~lhs, x, A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the default assignment kernel for the transpose sparse vector-
   // dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<VT1,VT2,MT1> >::Type
      selectAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      for( size_t j=0UL; j<N; ++j ) {
         y[j] = element->value() * A(element->index(),j);
      }

      ++element;

      for( ; element!=end; ++element ) {
         for( size_t j=0UL; j<N; ++j ) {
            y[j] += element->value() * A(element->index(),j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the optimized assignment kernel for the transpose sparse vector-
   // dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<VT1,VT2,MT1> >::Type
      selectAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      const size_t iend( x.nonZeros() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( x.nonZeros() - ( x.nonZeros() % 4UL ) ) == iend, "Invalid end calculation" );

      if( iend > 3UL )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );
         ++element;
         const size_t i2( element->index() );
         const VET    v2( element->value() );
         ++element;
         const size_t i3( element->index() );
         const VET    v3( element->value() );
         ++element;
         const size_t i4( element->index() );
         const VET    v4( element->value() );
         ++element;

         for( size_t j=0UL; j<N; ++j ) {
            y[j] = v1 * A(i1,j) + v2 * A(i2,j) + v3 * A(i3,j) + v4 * A(i4,j);
         }
      }
      else
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );
         ++element;

         for( size_t j=0UL; j<N; ++j ) {
            y[j] = v1 * A(i1,j);
         }
      }

      for( size_t i=(iend>3UL)?(4UL):(1UL); (i+4UL)<=iend; i+=4UL )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );
         ++element;
         const size_t i2( element->index() );
         const VET    v2( element->value() );
         ++element;
         const size_t i3( element->index() );
         const VET    v3( element->value() );
         ++element;
         const size_t i4( element->index() );
         const VET    v4( element->value() );
         ++element;

         for( size_t j=0UL; j<N; ++j ) {
            y[j] += v1 * A(i1,j) + v2 * A(i2,j) + v3 * A(i3,j) + v4 * A(i4,j);
         }
      }
      for( ; element!=end; ++element )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );

         for( size_t j=0UL; j<N; ++j ) {
            y[j] += v1 * A(i1,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized assignment to dense vectors******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the vectorized assignment kernel for the transpose sparse vector-
   // dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseVectorizedKernel<VT1,VT2,MT1> >::Type
      selectAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef IntrinsicTrait<ElementType>  IT;
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      const size_t iend( x.nonZeros() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( x.nonZeros() - ( x.nonZeros() % 4UL ) ) == iend, "Invalid end calculation" );

      if( iend > 3UL )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );
         ++element;
         const size_t        i2( element->index() );
         const IntrinsicType v2( set( element->value() ) );
         ++element;
         const size_t        i3( element->index() );
         const IntrinsicType v3( set( element->value() ) );
         ++element;
         const size_t        i4( element->index() );
         const IntrinsicType v4( set( element->value() ) );
         ++element;

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, v1 * A.load(i1,j) + v2 * A.load(i2,j) + v3 * A.load(i3,j) + v4 * A.load(i4,j) );
         }
      }
      else
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );
         ++element;

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, v1 * A.load(i1,j) );
         }
      }

      for( size_t i=(iend>3UL)?(4UL):(1UL); (i+4UL)<=iend; i+=4UL )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );
         ++element;
         const size_t        i2( element->index() );
         const IntrinsicType v2( set( element->value() ) );
         ++element;
         const size_t        i3( element->index() );
         const IntrinsicType v3( set( element->value() ) );
         ++element;
         const size_t        i4( element->index() );
         const IntrinsicType v4( set( element->value() ) );
         ++element;

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, y.load(j) + v1 * A.load(i1,j) + v2 * A.load(i2,j) + v3 * A.load(i3,j) + v4 * A.load(i4,j) );
         }
      }
      for( ; element!=end; ++element )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, y.load(j) + v1 * A.load(i1,j) );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the SMP assignment kernel for the transpose sparse vector-dense
   // matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,VT2,MT1> >::Type
      selectAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      smpAssign( y, x * A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse vector-dense matrix multiplication to a sparse
   //        vector (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse vector-
   // dense matrix multiplication expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT2,true>& lhs, const TSVecDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE  ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a transpose sparse vector-dense matrix multiplication to a
   //        dense vector (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose sparse
   // vector-dense matrix multiplication expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT2,true>& lhs, const TSVecDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the right-hand side dense matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse vector-dense matrix multiplication
      TSVecDMatMultExpr::selectAddAssignKernel( ~lhs, x, A );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the default addition assignment kernel for the transpose sparse
   // vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<VT1,VT2,MT1> >::Type
      selectAddAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      for( ; element!=end; ++element ) {
         for( size_t j=0UL; j<N; ++j ) {
            y[j] += element->value() * A(element->index(),j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the optimized addition assignment kernel for the transpose sparse
   // vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<VT1,VT2,MT1> >::Type
      selectAddAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      const size_t iend( x.nonZeros() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( x.nonZeros() - ( x.nonZeros() % 4UL ) ) == iend, "Invalid end calculation" );

      for( size_t i=0UL; (i+4UL)<=iend; i+=4UL )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );
         ++element;
         const size_t i2( element->index() );
         const VET    v2( element->value() );
         ++element;
         const size_t i3( element->index() );
         const VET    v3( element->value() );
         ++element;
         const size_t i4( element->index() );
         const VET    v4( element->value() );
         ++element;

         for( size_t j=0UL; j<N; ++j ) {
            y[j] += v1 * A(i1,j) + v2 * A(i2,j) + v3 * A(i3,j) + v4 * A(i4,j);
         }
      }
      for( ; element!=end; ++element )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );

         for( size_t j=0UL; j<N; ++j ) {
            y[j] += v1 * A(i1,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized addition assignment to dense vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized addition assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the vectorized addition assignment kernel for the transpose sparse
   // vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseVectorizedKernel<VT1,VT2,MT1> >::Type
      selectAddAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef IntrinsicTrait<ElementType>  IT;
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      const size_t iend( x.nonZeros() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( x.nonZeros() - ( x.nonZeros() % 4UL ) ) == iend, "Invalid end calculation" );

      for( size_t i=0UL; (i+4UL)<=iend; i+=4UL )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );
         ++element;
         const size_t        i2( element->index() );
         const IntrinsicType v2( set( element->value() ) );
         ++element;
         const size_t        i3( element->index() );
         const IntrinsicType v3( set( element->value() ) );
         ++element;
         const size_t        i4( element->index() );
         const IntrinsicType v4( set( element->value() ) );
         ++element;

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, y.load(j) + v1 * A.load(i1,j) + v2 * A.load(i2,j) + v3 * A.load(i3,j) + v4 * A.load(i4,j) );
         }
      }
      for( ; element!=end; ++element )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, y.load(j) + v1 * A.load(i1,j) );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the SMP addition assignment kernel for the transpose sparse vector-
   // dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,VT2,MT1> >::Type
      selectAddAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      smpAddAssign( y, x * A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a transpose sparse vector-dense matrix multiplication to a
   //        dense vector (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // sparse vector-dense matrix multiplication expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT2,true>& lhs, const TSVecDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the right-hand side dense matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse vector-dense matrix multiplication
      TSVecDMatMultExpr::selectSubAssignKernel( ~lhs, x, A );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose sparse vector-dense matrix
   //        multiplication (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the transpose
   // sparse vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<VT1,VT2,MT1> >::Type
      selectSubAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      for( ; element!=end; ++element ) {
         for( size_t j=0UL; j<N; ++j ) {
            y[j] -= element->value() * A(element->index(),j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to dense vectors*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a transpose sparse vector-dense matrix
   //        multiplication (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the optimized subtraction assignment kernel for the transpose
   // sparse vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<VT1,VT2,MT1> >::Type
      selectSubAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      const size_t iend( x.nonZeros() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( x.nonZeros() - ( x.nonZeros() % 4UL ) ) == iend, "Invalid end calculation" );

      for( size_t i=0UL; (i+4UL)<=iend; i+=4UL )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );
         ++element;
         const size_t i2( element->index() );
         const VET    v2( element->value() );
         ++element;
         const size_t i3( element->index() );
         const VET    v3( element->value() );
         ++element;
         const size_t i4( element->index() );
         const VET    v4( element->value() );
         ++element;

         for( size_t j=0UL; j<N; ++j ) {
            y[j] -= v1 * A(i1,j) + v2 * A(i2,j) + v3 * A(i3,j) + v4 * A(i4,j);
         }
      }
      for( ; element!=end; ++element )
      {
         const size_t i1( element->index() );
         const VET    v1( element->value() );

         for( size_t j=0UL; j<N; ++j ) {
            y[j] -= v1 * A(i1,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized subtraction assignment to dense vectors******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized subtraction assignment of a transpose sparse vector-dense matrix
   //        multiplication (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the vectorized subtraction assignment kernel for the transpose
   // sparse vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseVectorizedKernel<VT1,VT2,MT1> >::Type
      selectSubAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef IntrinsicTrait<ElementType>  IT;
      typedef typename RemoveReference<LT>::Type::ConstIterator  ConstIterator;

      BLAZE_INTERNAL_ASSERT( x.nonZeros() != 0UL, "Invalid number of non-zero elements" );

      const size_t N( A.columns() );

      ConstIterator element( x.begin() );
      const ConstIterator end( x.end() );

      const size_t iend( x.nonZeros() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( x.nonZeros() - ( x.nonZeros() % 4UL ) ) == iend, "Invalid end calculation" );

      for( size_t i=0UL; (i+4UL)<=iend; i+=4UL )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );
         ++element;
         const size_t        i2( element->index() );
         const IntrinsicType v2( set( element->value() ) );
         ++element;
         const size_t        i3( element->index() );
         const IntrinsicType v3( set( element->value() ) );
         ++element;
         const size_t        i4( element->index() );
         const IntrinsicType v4( set( element->value() ) );
         ++element;

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, y.load(j) - v1 * A.load(i1,j) - v2 * A.load(i2,j) - v3 * A.load(i3,j) - v4 * A.load(i4,j) );
         }
      }
      for( ; element!=x.end(); ++element )
      {
         const size_t        i1( element->index() );
         const IntrinsicType v1( set( element->value() ) );

         for( size_t j=0UL; j<N; j+=IT::size ) {
            y.store( j, y.load(j) - v1 * A.load(i1,j) );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose sparse vector-dense matrix multiplication
   //        (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the SMP subtraction assignment kernel for the transpose sparse
   // vector-dense matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,VT2,MT1> >::Type
      selectSubAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      smpSubAssign( y, x * A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a transpose sparse vector-dense matrix multiplication
   //        to a dense vector (\f$ \vec{y}^T*=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // sparse vector-dense matrix multiplication expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT2,true>& lhs, const TSVecDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE  ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
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
/*!\brief Multiplication operator for the multiplication of a transpose sparse vector and a
//        row-major dense matrix (\f$ \vec{y}^T=\vec{x}^T*A \f$).
// \ingroup dense_matrix
//
// \param vec The left-hand side transpose sparse vector for the multiplication.
// \param mat The right-hand side row-major dense matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose sparse vector and a row-major
// dense matrix:

   \code
   using blaze::rowVector;
   using blaze::rowMajor;

   blaze::CompressedVector<double,rowVector> x, y;
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// The operator returns an expression representing a transpose sparse vector of the higher-order
// element type of the two involved element types \a T1::ElementType and \a T2::ElementType.
// Both the dense matrix type \a T1 and the dense vector type \a T2 as well as the two element
// types \a T1::ElementType and \a T2::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of rows of
// the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1, typename T2 >
inline const typename DisableIf< IsMatMatMultExpr<T2>, TSVecDMatMultExpr<T1,T2> >::Type
   operator*( const SparseVector<T1,true>& vec, const DenseMatrix<T2,false>& mat )
{
   BLAZE_FUNCTION_TRACE;

   if( (~vec).size() != (~mat).rows() )
      throw std::invalid_argument( "Vector and matrix sizes do not match" );

   return TSVecDMatMultExpr<T1,T2>( ~vec, ~mat );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a transpose sparse vector and a
//        dense matrix-matrix multiplication expression (\f$ \vec{y}^T=\vec{x}^T*(A*B) \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side sparse vector for the multiplication.
// \param mat The right-hand side dense matrix-matrix multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a sparse
// vector and a dense matrix-matrix multiplication expression. It restructures the expression
// \f$ \vec{y}^T=\vec{x}^T*(A*B) \f$ to the expression \f$ \vec{y}^T=(\vec{x}^T*A)*B \f$.
*/
template< typename T1  // Type of the left-hand side sparse vector
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order of the right-hand side dense matrix
inline const typename EnableIf< IsMatMatMultExpr<T2>, MultExprTrait<T1,T2> >::Type::Type
   operator*( const SparseVector<T1,true>& vec, const DenseMatrix<T2,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( vec * (~mat).leftOperand() ) * (~mat).rightOperand();
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, bool AF >
struct SubvectorExprTrait< TSVecDMatMultExpr<VT,MT>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename MultExprTrait< VT, typename SubmatrixExprTrait<const MT,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
