//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatDMatMultExpr.h
//  \brief Header file for the sparse matrix/dense matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATDMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATDMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/DMatDVecMultExprTrait.h>
#include <blaze/math/traits/DMatSVecMultExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SMatDVecMultExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/TDVecDMatMultExprTrait.h>
#include <blaze/math/traits/TDVecSMatMultExprTrait.h>
#include <blaze/math/traits/TSVecDMatMultExprTrait.h>
#include <blaze/math/traits/TSVecSMatMultExprTrait.h>
#include <blaze/math/typetraits/Columns.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATDMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-dense matrix multiplications.
// \ingroup dense_matrix_expression
//
// The SMatDMatMultExpr class represents the compile time expression for multiplications between
// a row-major sparse matrix and a row-major dense matrix.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
class SMatDMatMultExpr : public DenseMatrix< SMatDMatMultExpr<MT1,MT2>, false >
                       , private MatMatMultExpr
                       , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT1::ResultType     RT1;  //!< Result type of the left-hand side sparse matrix expression.
   typedef typename MT2::ResultType     RT2;  //!< Result type of the right-hand side dense matrix expression.
   typedef typename RT1::ElementType    ET1;  //!< Element type of the left-hand side sparse matrix expression.
   typedef typename RT2::ElementType    ET2;  //!< Element type of the right-hand side dense matrix expression.
   typedef typename MT1::CompositeType  CT1;  //!< Composite type of the left-hand side sparse matrix expression.
   typedef typename MT2::CompositeType  CT2;  //!< Composite type of the right-hand side dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side sparse matrix expression.
   enum { evaluateLeft = IsComputation<MT1>::value || RequiresEvaluation<MT1>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   enum { evaluateRight = IsComputation<MT2>::value || RequiresEvaluation<MT2>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The CanExploitSymmetry struct is a helper struct for the selection of the optimal
       evaluation strategy. In case the target matrix is column-major and either of the
       two matrix operands is symmetric, \a value is set to 1 and an optimized evaluation
       strategy is selected. Otherwise \a value is set to 0 and the default strategy is
       chosen. */
   template< typename T1, typename T2, typename T3 >
   struct CanExploitSymmetry {
      enum { value = IsColumnMajorMatrix<T1>::value &&
                     ( IsSymmetric<T2>::value || IsSymmetric<T3>::value ) };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The IsEvaluationRequired struct is a helper struct for the selection of the parallel
       evaluation strategy. In case either of the two matrix operands requires an intermediate
       evaluation, the nested \value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct IsEvaluationRequired {
      enum { value = ( evaluateLeft || evaluateRight ) &&
                     !CanExploitSymmetry<T1,T2,T3>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case all three involved data types are suited for a vectorized computation of the
       matrix multiplication, the nested \value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseVectorizedKernel {
      enum { value = !IsDiagonal<T3>::value &&
                     T1::vectorizable && T3::vectorizable &&
                     IsRowMajorMatrix<T1>::value &&
                     IsSame<typename T1::ElementType,typename T2::ElementType>::value &&
                     IsSame<typename T1::ElementType,typename T3::ElementType>::value &&
                     IntrinsicTrait<typename T1::ElementType>::addition &&
                     IntrinsicTrait<typename T1::ElementType>::subtraction &&
                     IntrinsicTrait<typename T1::ElementType>::multiplication };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case a vectorized computation of the matrix multiplication is not possible, but a
       loop-unrolled computation is feasible, the nested \value will be set to 1, otherwise
       it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseOptimizedKernel {
      enum { value = !UseVectorizedKernel<T1,T2,T3>::value &&
                     !IsDiagonal<T3>::value &&
                     !IsResizable<typename T1::ElementType>::value &&
                     !IsResizable<ET1>::value };
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
      enum { value = !UseVectorizedKernel<T1,T2,T3>::value &&
                     !UseOptimizedKernel<T1,T2,T3>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SMatDMatMultExpr<MT1,MT2>                   This;           //!< Type of this SMatDMatMultExpr instance.
   typedef typename MultTrait<RT1,RT2>::Type           ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType           OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const ElementType                           ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                            CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT1>::value, const MT1, const MT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT2>::value, const MT2, const MT2& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side sparse matrix operand.
   typedef typename SelectType< evaluateLeft, const RT1, CT1 >::Type  LT;

   //! Type for the assignment of the right-hand side dense matrix operand.
   typedef typename SelectType< evaluateRight, const RT2, CT2 >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = !IsDiagonal<MT2>::value &&
                         MT2::vectorizable &&
                         IsSame<ET1,ET2>::value &&
                         IntrinsicTrait<ET1>::addition &&
                         IntrinsicTrait<ET1>::multiplication };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateLeft  && MT1::smpAssignable &&
                          !evaluateRight && MT2::smpAssignable };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatDMatMultExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the multiplication expression.
   // \param rhs The right-hand side dense matrix operand of the multiplication expression.
   */
   explicit inline SMatDMatMultExpr( const MT1& lhs, const MT2& rhs )
      : lhs_( lhs )  // Left-hand side sparse matrix of the multiplication expression
      , rhs_( rhs )  // Right-hand side dense matrix of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.columns() == rhs.rows(), "Invalid matrix sizes" );
   }
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < lhs_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < rhs_.columns(), "Invalid column access index" );

      typedef typename RemoveReference<CT1>::Type::ConstIterator  ConstIterator;

      ElementType tmp = ElementType();

      // Early exit
      if( lhs_.columns() == 0UL )
         return tmp;

      // Fast computation in case the left-hand side sparse matrix directly provides iterators
      if( !RequiresEvaluation<MT1>::value )
      {
         CT1 A( lhs_ );  // Evaluation of the left-hand side sparse matrix operand

         const ConstIterator end( ( IsUpper<MT2>::value )
                                  ?( IsStrictlyUpper<MT2>::value ? A.lowerBound(i,j) : A.upperBound(i,j) )
                                  :( A.end(i) ) );
         ConstIterator element( ( IsLower<MT2>::value )
                                ?( IsStrictlyLower<MT2>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                :( A.begin(i) ) );

         if( element != end ) {
            tmp = element->value() * rhs_(element->index(),j);
            ++element;
            for( ; element!=end; ++element ) {
               tmp += element->value() * rhs_(element->index(),j);
            }
         }
      }

      // Default computation in case the left-hand side sparse matrix doesn't provide iterators
      else
      {
         const size_t kbegin( ( IsUpper<MT1>::value )
                              ?( ( IsLower<MT2>::value )
                                 ?( max( ( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                                       , ( IsStrictlyLower<MT2>::value ? j+1UL : j ) ) )
                                 :( IsStrictlyUpper<MT1>::value ? i+1UL : i ) )
                              :( ( IsLower<MT2>::value )
                                 ?( IsStrictlyLower<MT2>::value ? j+1UL : j )
                                 :( 0UL ) ) );
         const size_t kend( ( IsLower<MT1>::value )
                            ?( ( IsUpper<MT2>::value )
                               ?( min( ( IsStrictlyLower<MT1>::value ? i : i+1UL )
                                     , ( IsStrictlyUpper<MT2>::value ? j : j+1UL ) ) )
                               :( IsStrictlyLower<MT1>::value ? i : i+1UL ) )
                            :( ( IsUpper<MT2>::value )
                               ?( IsStrictlyUpper<MT2>::value ? j : j+1UL )
                               :( lhs_.columns() ) ) );

         if( ( !IsTriangular<MT1>::value && !IsTriangular<MT2>::value ) || kbegin < kend ) {
            tmp = lhs_(i,kbegin) * rhs_(kbegin,j);
            for( size_t k=kbegin+1UL; k<kend; ++k ) {
               tmp += lhs_(i,k) * rhs_(k,j);
            }
         }
      }

      return tmp;
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const {
      return lhs_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const {
      return rhs_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse matrix operand.
   //
   // \return The left-hand side sparse matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense matrix operand.
   //
   // \return The right-hand side dense matrix operand.
   */
   inline RightOperand rightOperand() const {
      return rhs_;
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
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
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

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const {
      return rhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return ( rows() > SMP_SMATDMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-dense matrix multiplication to a dense matrix
   //        (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-dense
   // matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline typename DisableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      assign( DenseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      SMatDMatMultExpr::selectAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to row-major dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a sparse matrix-dense matrix multiplication to row-major
   //        dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the default row-major assignment kernel for the sparse matrix-
   // dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<MT3,MT4,MT5> >::Type
      selectAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         for( size_t j=0UL; j<(~C).columns(); ++j ) {
            reset( (~C)(i,j) );
         }

         ConstIterator element( A.begin(i) );
         const ConstIterator end( A.end(i) );

         for( ; element!=end; ++element )
         {
            const size_t index( element->index() );

            if( IsDiagonal<MT5>::value )
            {
               (~C)(i,index) = element->value() * B(index,index);
            }
            else
            {
               const size_t jbegin( ( IsUpper<MT5>::value )
                                    ?( IsStrictlyUpper<MT5>::value ? index+1UL : index )
                                    :( 0UL ) );
               const size_t jend( ( IsLower<MT5>::value )
                                  ?( IsStrictlyLower<MT5>::value ? index : index+1UL )
                                  :( B.columns() ) );
               BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

               for( size_t j=jbegin; j<jend; ++j ) {
                  if( isDefault( (~C)(i,j) ) )
                     (~C)(i,j) = element->value() * B(index,j);
                  else
                     (~C)(i,j) += element->value() * B(index,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to row-major dense matrices********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a sparse matrix-dense matrix multiplication to row-major
   //        dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the optimized row-major assignment kernel for the sparse matrix-
   // dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<MT3,MT4,MT5> >::Type
      selectAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      reset( ~C );

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         const size_t nonzeros( A.nonZeros(i) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t i1( element->index() );
            const ET1    v1( element->value() );
            ++element;
            const size_t i2( element->index() );
            const ET1    v2( element->value() );
            ++element;
            const size_t i3( element->index() );
            const ET1    v3( element->value() );
            ++element;
            const size_t i4( element->index() );
            const ET1    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( i1 < i2 && i2 < i3 && i3 < i4, "Invalid sparse matrix index detected" );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i4 : i4+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            const size_t jnum( jend - jbegin );
            const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

            for( size_t j=jbegin; j<jpos; j+=4UL ) {
               (~C)(i,j    ) += v1 * B(i1,j    ) + v2 * B(i2,j    ) + v3 * B(i3,j    ) + v4 * B(i4,j    );
               (~C)(i,j+1UL) += v1 * B(i1,j+1UL) + v2 * B(i2,j+1UL) + v3 * B(i3,j+1UL) + v4 * B(i4,j+1UL);
               (~C)(i,j+2UL) += v1 * B(i1,j+2UL) + v2 * B(i2,j+2UL) + v3 * B(i3,j+2UL) + v4 * B(i4,j+2UL);
               (~C)(i,j+3UL) += v1 * B(i1,j+3UL) + v2 * B(i2,j+3UL) + v3 * B(i3,j+3UL) + v4 * B(i4,j+3UL);
            }
            for( size_t j=jpos; j<jend; ++j ) {
               (~C)(i,j) += v1 * B(i1,j) + v2 * B(i2,j) + v3 * B(i3,j) + v4 * B(i4,j);
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t i1( element->index() );
            const ET1    v1( element->value() );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            const size_t jnum( jend - jbegin );
            const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

            for( size_t j=jbegin; j<jpos; j+=4UL ) {
               (~C)(i,j    ) += v1 * B(i1,j    );
               (~C)(i,j+1UL) += v1 * B(i1,j+1UL);
               (~C)(i,j+2UL) += v1 * B(i1,j+2UL);
               (~C)(i,j+3UL) += v1 * B(i1,j+3UL);
            }
            for( size_t j=jpos; j<jend; ++j ) {
               (~C)(i,j) += v1 * B(i1,j);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized assignment to row-major dense matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized assignment of a sparse matrix-dense matrix multiplication to row-major
   //        dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the vectorized row-major assignment kernel for the sparse matrix-
   // dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseVectorizedKernel<MT3,MT4,MT5> >::Type
      selectAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef IntrinsicTrait<ElementType>  IT;
      typedef typename MT4::ConstIterator  ConstIterator;

      reset( ~C );

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         const size_t nonzeros( A.nonZeros(i) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
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

            BLAZE_INTERNAL_ASSERT( i1 < i2 && i2 < i3 && i3 < i4, "Invalid sparse matrix index detected" );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( ( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ) & size_t(-IT::size) )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i4 : i4+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; j+=IT::size ) {
               (~C).store( i, j, (~C).load(i,j) + v1 * B.load(i1,j) + v2 * B.load(i2,j) + v3 * B.load(i3,j) + v4 * B.load(i4,j) );
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t        i1( element->index() );
            const IntrinsicType v1( set( element->value() ) );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( ( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ) & size_t(-IT::size) )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; j+=IT::size ) {
               (~C).store( i, j, (~C).load(i,j) + v1 * B.load(i1,j) );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to column-major dense matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a sparse matrix-dense matrix multiplication to column-major
   //        dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default column-major assignment kernel for the sparse matrix-
   // dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAssignKernel( DenseMatrix<MT3,true>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      const size_t jpos( B.columns() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( B.columns() - ( B.columns() % 4UL ) ) == jpos, "Invalid end calculation" );

      for( size_t j=0UL; j<jpos; j+=4UL ) {
         for( size_t i=0UL; i<A.rows(); ++i )
         {
            ConstIterator element( ( IsLower<MT5>::value )
                                   ?( IsStrictlyLower<MT5>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                   :( A.begin(i) ) );
            const ConstIterator end( ( IsUpper<MT5>::value )
                                     ?( IsStrictlyUpper<MT5>::value ? A.lowerBound(i,j+4UL) : A.upperBound(i,j+4UL) )
                                     :( A.end(i) ) );

            if( element == end ) {
               reset( (~C)(i,j    ) );
               reset( (~C)(i,j+1UL) );
               reset( (~C)(i,j+2UL) );
               reset( (~C)(i,j+3UL) );
               continue;
            }

            (~C)(i,j    ) = element->value() * B(element->index(),j    );
            (~C)(i,j+1UL) = element->value() * B(element->index(),j+1UL);
            (~C)(i,j+2UL) = element->value() * B(element->index(),j+2UL);
            (~C)(i,j+3UL) = element->value() * B(element->index(),j+3UL);
            ++element;
            for( ; element!=end; ++element ) {
               (~C)(i,j    ) += element->value() * B(element->index(),j    );
               (~C)(i,j+1UL) += element->value() * B(element->index(),j+1UL);
               (~C)(i,j+2UL) += element->value() * B(element->index(),j+2UL);
               (~C)(i,j+3UL) += element->value() * B(element->index(),j+3UL);
            }
         }
      }

      for( size_t j=jpos; j<B.columns(); ++j ) {
         for( size_t i=0UL; i<A.rows(); ++i )
         {
            ConstIterator element( ( IsLower<MT5>::value )
                                   ?( IsStrictlyLower<MT5>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                   :( A.begin(i) ) );
            const ConstIterator end( ( IsUpper<MT5>::value )
                                     ?( IsStrictlyUpper<MT5>::value ? A.lowerBound(i,j) : A.upperBound(i,j) )
                                     :( A.end(i) ) );

            if( element == end ) {
               reset( (~C)(i,j) );
               continue;
            }

            (~C)(i,j) = element->value() * B(element->index(),j);
            ++element;
            for( ; element!=end; ++element )
               (~C)(i,j) += element->value() * B(element->index(),j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-dense matrix multiplication to a sparse matrix
   //        (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-dense
   // matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      typedef typename SelectType< SO, OppositeType, ResultType >::Type  TmpType;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename TmpType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( serial( rhs ) );
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring assignment to column-major matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring assignment of a sparse matrix-dense matrix multiplication to a
   //        column-major matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring assignment of a sparse matrix-
   // dense matrix multiplication expression to a column-major matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline typename EnableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      assign( Matrix<MT,true>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( IsSymmetric<MT1>::value )
         assign( ~lhs, trans( rhs.lhs_ ) * rhs.rhs_ );
      else
         assign( ~lhs, rhs.lhs_ * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-dense matrix multiplication to a dense matrix
   //        (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      SMatDMatMultExpr::selectAddAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to row-major dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a sparse matrix-dense matrix multiplication to
   //        row-major dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the default row-major addition assignment kernel for the sparse
   // matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<MT3,MT4,MT5> >::Type
      selectAddAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         for( ; element!=end; ++element )
         {
            const size_t i1( element->index() );

            if( IsDiagonal<MT5>::value )
            {
               (~C)(i,i1) += element->value() * B(i1,i1);
            }
            else
            {
               const size_t jbegin( ( IsUpper<MT5>::value )
                                    ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ):( 0UL ) );
               const size_t jend( ( IsLower<MT5>::value )
                                  ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                                  :( B.columns() ) );
               BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

               const size_t jnum( jend - jbegin );
               const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
               BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

               for( size_t j=jbegin; j<jpos; j+=4UL ) {
                  (~C)(i,j    ) += element->value() * B(i1,j    );
                  (~C)(i,j+1UL) += element->value() * B(i1,j+1UL);
                  (~C)(i,j+2UL) += element->value() * B(i1,j+2UL);
                  (~C)(i,j+3UL) += element->value() * B(i1,j+3UL);
               }
               for( size_t j=jpos; j<jend; ++j ) {
                  (~C)(i,j) += element->value() * B(i1,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to row-major dense matrices***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a sparse matrix-dense matrix multiplication to
   //        row-major dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the optimized row-major addition assignment kernel for the sparse
   // matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<MT3,MT4,MT5> >::Type
      selectAddAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         const size_t nonzeros( A.nonZeros(i) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t i1( element->index() );
            const ET1    v1( element->value() );
            ++element;
            const size_t i2( element->index() );
            const ET1    v2( element->value() );
            ++element;
            const size_t i3( element->index() );
            const ET1    v3( element->value() );
            ++element;
            const size_t i4( element->index() );
            const ET1    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( i1 < i2 && i2 < i3 && i3 < i4, "Invalid sparse matrix index detected" );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i4 : i4+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            const size_t jnum( jend - jbegin );
            const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

            for( size_t j=jbegin; j<jpos; j+=4UL ) {
               (~C)(i,j    ) += v1 * B(i1,j    ) + v2 * B(i2,j    ) + v3 * B(i3,j    ) + v4 * B(i4,j    );
               (~C)(i,j+1UL) += v1 * B(i1,j+1UL) + v2 * B(i2,j+1UL) + v3 * B(i3,j+1UL) + v4 * B(i4,j+1UL);
               (~C)(i,j+2UL) += v1 * B(i1,j+2UL) + v2 * B(i2,j+2UL) + v3 * B(i3,j+2UL) + v4 * B(i4,j+2UL);
               (~C)(i,j+3UL) += v1 * B(i1,j+3UL) + v2 * B(i2,j+3UL) + v3 * B(i3,j+3UL) + v4 * B(i4,j+3UL);
            }
            for( size_t j=jpos; j<jend; ++j ) {
               (~C)(i,j) += v1 * B(i1,j) + v2 * B(i2,j) + v3 * B(i3,j) + v4 * B(i4,j);
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t i1( element->index() );
            const ET1    v1( element->value() );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            const size_t jnum( jend - jbegin );
            const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

            for( size_t j=jbegin; j<jpos; j+=4UL ) {
               (~C)(i,j    ) += v1 * B(i1,j    );
               (~C)(i,j+1UL) += v1 * B(i1,j+1UL);
               (~C)(i,j+2UL) += v1 * B(i1,j+2UL);
               (~C)(i,j+3UL) += v1 * B(i1,j+3UL);
            }
            for( size_t j=jpos; j<jend; ++j ) {
               (~C)(i,j) += v1 * B(i1,j);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized addition assignment to row-major dense matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized addition assignment of a sparse matrix-dense matrix multiplication to
   //        row-major dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the vectorized row-major addition assignment kernel for the sparse
   // matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseVectorizedKernel<MT3,MT4,MT5> >::Type
      selectAddAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef IntrinsicTrait<ElementType>  IT;
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         const size_t nonzeros( A.nonZeros(i) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
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

            BLAZE_INTERNAL_ASSERT( i1 < i2 && i2 < i3 && i3 < i4, "Invalid sparse matrix index detected" );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( ( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ) & size_t(-IT::size) )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i4 : i4+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; j+=IT::size ) {
               (~C).store( i, j, (~C).load(i,j) + v1 * B.load(i1,j) + v2 * B.load(i2,j) + v3 * B.load(i3,j) + v4 * B.load(i4,j) );
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t        i1( element->index() );
            const IntrinsicType v1( set( element->value() ) );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( ( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ) & size_t(-IT::size) )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; j+=IT::size ) {
               (~C).store( i, j, (~C).load(i,j) + v1 * B.load(i1,j) );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to column-major dense matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a sparse matrix-dense matrix multiplication to
   //        column-major dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default column-major addition assignment kernel for the sparse
   // matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAddAssignKernel( DenseMatrix<MT3,true>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      const size_t jpos( B.columns() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( B.columns() - ( B.columns() % 4UL ) ) == jpos, "Invalid end calculation" );

      for( size_t j=0UL; j<jpos; j+=4UL ) {
         for( size_t i=0UL; i<A.rows(); ++i )
         {
            ConstIterator element( ( IsLower<MT5>::value )
                                   ?( IsStrictlyLower<MT5>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                   :( A.begin(i) ) );
            const ConstIterator end( ( IsUpper<MT5>::value )
                                     ?( IsStrictlyUpper<MT5>::value ? A.lowerBound(i,j+4UL) : A.upperBound(i,j+4UL) )
                                     :( A.end(i) ) );

            for( ; element!=end; ++element ) {
               (~C)(i,j    ) += element->value() * B(element->index(),j    );
               (~C)(i,j+1UL) += element->value() * B(element->index(),j+1UL);
               (~C)(i,j+2UL) += element->value() * B(element->index(),j+2UL);
               (~C)(i,j+3UL) += element->value() * B(element->index(),j+3UL);
            }
         }
      }

      for( size_t j=jpos; j<B.columns(); ++j ) {
         for( size_t i=0UL; i<A.rows(); ++i )
         {
            ConstIterator element( ( IsLower<MT5>::value )
                                   ?( IsStrictlyLower<MT5>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                   :( A.begin(i) ) );
            const ConstIterator end( ( IsUpper<MT5>::value )
                                     ?( IsStrictlyUpper<MT5>::value ? A.lowerBound(i,j) : A.upperBound(i,j) )
                                     :( A.end(i) ) );

            for( ; element!=end; ++element )
               (~C)(i,j) += element->value() * B(element->index(),j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring addition assignment to column-major matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring addition assignment of a sparse matrix-dense matrix multiplication to
   //        a column-major matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring addition assignment of a sparse
   // matrix-dense matrix multiplication expression to a column-major matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler in
   // case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline typename EnableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      addAssign( Matrix<MT,true>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( IsSymmetric<MT1>::value )
         addAssign( ~lhs, trans( rhs.lhs_ ) * rhs.rhs_ );
      else
         addAssign( ~lhs, rhs.lhs_ * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-dense matrix multiplication to a dense
   //        matrix (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      SMatDMatMultExpr::selectSubAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to row-major dense matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a sparse matrix-dense matrix multiplication
   //        (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the default row-major subtraction assignment kernel for the sparse
   // matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<MT3,MT4,MT5> >::Type
      selectSubAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         for( ; element!=end; ++element )
         {
            const size_t i1( element->index() );

            if( IsDiagonal<MT5>::value )
            {
               (~C)(i,i1) -= element->value() * B(i1,i1);
            }
            else
            {
               const size_t jbegin( ( IsUpper<MT5>::value )
                                    ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ):( 0UL ) );
               const size_t jend( ( IsLower<MT5>::value )
                                  ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                                  :( B.columns() ) );
               BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

               const size_t jnum( jend - jbegin );
               const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
               BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

               for( size_t j=jbegin; j<jpos; j+=4UL ) {
                  (~C)(i,j    ) -= element->value() * B(i1,j    );
                  (~C)(i,j+1UL) -= element->value() * B(i1,j+1UL);
                  (~C)(i,j+2UL) -= element->value() * B(i1,j+2UL);
                  (~C)(i,j+3UL) -= element->value() * B(i1,j+3UL);
               }
               for( size_t j=jpos; j<jend; ++j ) {
                  (~C)(i,j) -= element->value() * B(i1,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to row-major dense matrices********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a sparse matrix-dense matrix multiplication
   //        (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the optimized row-major subtraction assignment kernel for the
   // sparse matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<MT3,MT4,MT5> >::Type
      selectSubAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         const size_t nonzeros( A.nonZeros(i) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t i1( element->index() );
            const ET1    v1( element->value() );
            ++element;
            const size_t i2( element->index() );
            const ET1    v2( element->value() );
            ++element;
            const size_t i3( element->index() );
            const ET1    v3( element->value() );
            ++element;
            const size_t i4( element->index() );
            const ET1    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( i1 < i2 && i2 < i3 && i3 < i4, "Invalid sparse matrix index detected" );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i4 : i4+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            const size_t jnum( jend - jbegin );
            const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

            for( size_t j=jbegin; j<jpos; j+=4UL ) {
               (~C)(i,j    ) -= v1 * B(i1,j    ) + v2 * B(i2,j    ) + v3 * B(i3,j    ) + v4 * B(i4,j    );
               (~C)(i,j+1UL) -= v1 * B(i1,j+1UL) + v2 * B(i2,j+1UL) + v3 * B(i3,j+1UL) + v4 * B(i4,j+1UL);
               (~C)(i,j+2UL) -= v1 * B(i1,j+2UL) + v2 * B(i2,j+2UL) + v3 * B(i3,j+2UL) + v4 * B(i4,j+2UL);
               (~C)(i,j+3UL) -= v1 * B(i1,j+3UL) + v2 * B(i2,j+3UL) + v3 * B(i3,j+3UL) + v4 * B(i4,j+3UL);
            }
            for( size_t j=jpos; j<jend; ++j ) {
               (~C)(i,j) -= v1 * B(i1,j) + v2 * B(i2,j) + v3 * B(i3,j) + v4 * B(i4,j);
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t i1( element->index() );
            const ET1    v1( element->value() );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            const size_t jnum( jend - jbegin );
            const size_t jpos( jbegin + ( jnum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( jbegin + jnum - ( jnum % 4UL ) ) == jpos, "Invalid end calculation" );

            for( size_t j=jbegin; j<jpos; j+=4UL ) {
               (~C)(i,j    ) -= v1 * B(i1,j    );
               (~C)(i,j+1UL) -= v1 * B(i1,j+1UL);
               (~C)(i,j+2UL) -= v1 * B(i1,j+2UL);
               (~C)(i,j+3UL) -= v1 * B(i1,j+3UL);
            }
            for( size_t j=jpos; j<jend; ++j ) {
               (~C)(i,j) -= v1 * B(i1,j);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized subtraction assignment to row-major dense matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized subtraction assignment of a sparse matrix-dense matrix multiplication
   //        (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side dense matrix operand.
   // \return void
   //
   // This function implements the vectorized row-major subtraction assignment kernel for the
   // sparse matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseVectorizedKernel<MT3,MT4,MT5> >::Type
      selectSubAssignKernel( DenseMatrix<MT3,false>& C, const MT4& A, const MT5& B )
   {
      typedef IntrinsicTrait<ElementType>  IT;
      typedef typename MT4::ConstIterator  ConstIterator;

      for( size_t i=0UL; i<A.rows(); ++i )
      {
         const ConstIterator end( A.end(i) );
         ConstIterator element( A.begin(i) );

         const size_t nonzeros( A.nonZeros(i) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
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

            BLAZE_INTERNAL_ASSERT( i1 < i2 && i2 < i3 && i3 < i4, "Invalid sparse matrix index detected" );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( ( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ) & size_t(-IT::size) )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i4 : i4+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; j+=IT::size ) {
               (~C).store( i, j, (~C).load(i,j) - v1 * B.load(i1,j) - v2 * B.load(i2,j) - v3 * B.load(i3,j) - v4 * B.load(i4,j) );
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t        i1( element->index() );
            const IntrinsicType v1( set( element->value() ) );

            const size_t jbegin( ( IsUpper<MT5>::value )
                                 ?( ( IsStrictlyUpper<MT5>::value ? i1+1UL : i1 ) & size_t(-IT::size) )
                                 :( 0UL ) );
            const size_t jend( ( IsLower<MT5>::value )
                               ?( IsStrictlyLower<MT5>::value ? i1 : i1+1UL )
                               :( B.columns() ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; j+=IT::size ) {
               (~C).store( i, j, (~C).load(i,j) - v1 * B.load(i1,j) );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to column-major dense matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a sparse matrix-dense matrix multiplication to
   //        column-major matrices (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default column-major subtraction assignment kernel for the
   // sparse matrix-dense matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectSubAssignKernel( DenseMatrix<MT3,true>& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      const size_t jpos( B.columns() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( B.columns() - ( B.columns() % 4UL ) ) == jpos, "Invalid end calculation" );

      for( size_t j=0UL; j<jpos; j+=4UL ) {
         for( size_t i=0UL; i<A.rows(); ++i )
         {
            ConstIterator element( ( IsLower<MT5>::value )
                                   ?( IsStrictlyLower<MT5>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                   :( A.begin(i) ) );
            const ConstIterator end( ( IsUpper<MT5>::value )
                                     ?( IsStrictlyUpper<MT5>::value ? A.lowerBound(i,j+4UL) : A.upperBound(i,j+4UL) )
                                     :( A.end(i) ) );

            for( ; element!=end; ++element ) {
               (~C)(i,j    ) -= element->value() * B(element->index(),j    );
               (~C)(i,j+1UL) -= element->value() * B(element->index(),j+1UL);
               (~C)(i,j+2UL) -= element->value() * B(element->index(),j+2UL);
               (~C)(i,j+3UL) -= element->value() * B(element->index(),j+3UL);
            }
         }
      }

      for( size_t j=jpos; j<B.columns(); ++j ) {
         for( size_t i=0UL; i<A.rows(); ++i )
         {
            ConstIterator element( ( IsLower<MT5>::value )
                                   ?( IsStrictlyLower<MT5>::value ? A.upperBound(i,j) : A.lowerBound(i,j) )
                                   :( A.begin(i) ) );
            const ConstIterator end( ( IsUpper<MT5>::value )
                                     ?( IsStrictlyUpper<MT5>::value ? A.lowerBound(i,j) : A.upperBound(i,j) )
                                     :( A.end(i) ) );

            for( ; element!=end; ++element )
               (~C)(i,j) -= element->value() * B(element->index(),j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring subtraction assignment to column-major matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring subtraction assignment of a sparse matrix-dense matrix multiplication
   //        to a column-major matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring subtraction assignment of a sparse
   // matrix-dense matrix multiplication expression to a column-major matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler in
   // case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline typename EnableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      subAssign( Matrix<MT,true>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( IsSymmetric<MT1>::value )
         subAssign( ~lhs, trans( rhs.lhs_ ) * rhs.rhs_ );
      else
         subAssign( ~lhs, rhs.lhs_ * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse matrix-dense matrix multiplication to a dense matrix
   //        (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix-dense
   // matrix multiplication expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle this function can only be selected by the compiler in case either of the
   // two matrix operands requires an intermediate evaluation and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline typename EnableIf< IsEvaluationRequired<MT,MT1,MT2> >::Type
      smpAssign( DenseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      smpAssign( ~lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse matrix-dense matrix multiplication to a sparse matrix
   //        (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix-dense
   // matrix multiplication expression to a sparse matrix. Due to the explicit application of the
   // SFINAE principle this function can only be selected by the compiler in case either of the
   // two matrix operands requires an intermediate evaluation and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline typename EnableIf< IsEvaluationRequired<MT,MT1,MT2> >::Type
      smpAssign( SparseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      typedef typename SelectType< SO, OppositeType, ResultType >::Type  TmpType;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename TmpType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( rhs );
      smpAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP assignment to column-major matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP assignment of a sparse matrix-dense matrix multiplication to a
   //        column-major matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP assignment of a sparse
   // matrix-dense matrix multiplication expression to a column-major matrix. Due to the
   // explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline typename EnableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      smpAssign( Matrix<MT,true>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( IsSymmetric<MT1>::value )
         smpAssign( ~lhs, trans( rhs.lhs_ ) * rhs.rhs_ );
      else
         smpAssign( ~lhs, rhs.lhs_ * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix-dense matrix multiplication to a dense
   //        matrix (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix-dense matrix multiplication expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case either of the two matrix operands requires an intermediate evaluation and no
   // symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline typename EnableIf< IsEvaluationRequired<MT,MT1,MT2> >::Type
      smpAddAssign( DenseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      smpAddAssign( ~lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP addition assignment to column-major matrices******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP addition assignment of a sparse matrix-dense matrix multiplication
   //        to a column-major matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP addition assignment of a sparse
   // matrix-dense matrix multiplication expression to a column-major matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler in
   // case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline typename EnableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      smpAddAssign( Matrix<MT,true>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( IsSymmetric<MT1>::value )
         smpAddAssign( ~lhs, trans( rhs.lhs_ ) * rhs.rhs_ );
      else
         smpAddAssign( ~lhs, rhs.lhs_ * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix-dense matrix multiplication to a dense
   //        matrix (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix-dense matrix multiplication expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case either of the two matrix operands requires an intermediate evaluation and no
   // symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline typename EnableIf< IsEvaluationRequired<MT,MT1,MT2> >::Type
      smpSubAssign( DenseMatrix<MT,SO>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      smpSubAssign( ~lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP subtraction assignment to column-major matrices***************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP subtraction assignment of a sparse matrix-dense matrix multiplication
   //        to a column-major matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP subtraction assignment of a
   // sparse matrix-dense matrix multiplication expression to a column-major matrix. Due to the
   // explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline typename EnableIf< CanExploitSymmetry<MT,MT1,MT2> >::Type
      smpSubAssign( Matrix<MT,true>& lhs, const SMatDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( IsSymmetric<MT1>::value )
         smpSubAssign( ~lhs, trans( rhs.lhs_ ) * rhs.rhs_ );
      else
         smpSubAssign( ~lhs, rhs.lhs_ * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATMULTEXPR( MT1, MT2 );
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
/*!\brief Multiplication operator for the multiplication of a row-major sparse matrix and a
//        row-major dense matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the multiplication.
// \param rhs The right-hand side dense matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a row-major sparse matrix and a row-major
// dense matrix:

   \code
   using blaze::rowMajor;

   blaze::CompressedMatrix<double,rowMajor> A;
   blaze::DynamicMatrix<double,rowMajor> B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a T1::ElementType and \a T2::ElementType.
// Both matrix types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side dense matrix
inline const SMatDMatMultExpr<T1,T2>
   operator*( const SparseMatrix<T1,false>& lhs, const DenseMatrix<T2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   return SMatDMatMultExpr<T1,T2>( ~lhs, ~rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  ROWS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct Rows< SMatDMatMultExpr<MT1,MT2> >
   : public Rows<MT1>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct Columns< SMatDMatMultExpr<MT1,MT2> >
   : public Columns<MT2>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsLower< SMatDMatMultExpr<MT1,MT2> >
   : public IsTrue< And< IsLower<MT1>, IsLower<MT2> >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsUniLower< SMatDMatMultExpr<MT1,MT2> >
   : public IsTrue< And< IsUniLower<MT1>, IsUniLower<MT2> >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsStrictlyLower< SMatDMatMultExpr<MT1,MT2> >
   : public IsTrue< Or< And< IsStrictlyLower<MT1>, IsLower<MT2> >
                      , And< IsStrictlyLower<MT2>, IsLower<MT1> > >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsUpper< SMatDMatMultExpr<MT1,MT2> >
   : public IsTrue< And< IsUpper<MT1>, IsUpper<MT2> >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsUniUpper< SMatDMatMultExpr<MT1,MT2> >
   : public IsTrue< And< IsUniUpper<MT1>, IsUniUpper<MT2> >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsStrictlyUpper< SMatDMatMultExpr<MT1,MT2> >
   : public IsTrue< Or< And< IsStrictlyUpper<MT1>, IsUpper<MT2> >
                      , And< IsStrictlyUpper<MT2>, IsUpper<MT1> > >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename VT >
struct DMatDVecMultExprTrait< SMatDMatMultExpr<MT1,MT2>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value &&
                                IsDenseVector<VT>::value   && IsColumnVector<VT>::value
                              , typename SMatDVecMultExprTrait< MT1, typename DMatDVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename VT >
struct DMatSVecMultExprTrait< SMatDMatMultExpr<MT1,MT2>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value &&
                                IsSparseVector<VT>::value  && IsColumnVector<VT>::value
                              , typename SMatDVecMultExprTrait< MT1, typename DMatSVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT1, typename MT2 >
struct TDVecDMatMultExprTrait< VT, SMatDMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value   && IsRowVector<VT>::value       &&
                                IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value
                              , typename TDVecDMatMultExprTrait< typename TDVecSMatMultExprTrait<VT,MT1>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT1, typename MT2 >
struct TSVecDMatMultExprTrait< VT, SMatDMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value  && IsRowVector<VT>::value       &&
                                IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value
                              , typename TSVecDMatMultExprTrait< typename TSVecSMatMultExprTrait<VT,MT1>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool AF >
struct SubmatrixExprTrait< SMatDMatMultExpr<MT1,MT2>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename MultExprTrait< typename SubmatrixExprTrait<const MT1,AF>::Type
                                 , typename SubmatrixExprTrait<const MT2,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct RowExprTrait< SMatDMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename MultExprTrait< typename RowExprTrait<const MT1>::Type, MT2 >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct ColumnExprTrait< SMatDMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename MultExprTrait< MT1, typename ColumnExprTrait<const MT2>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
