//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatDVecMultExpr.h
//  \brief Header file for the transpose dense matrix/dense vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/gemv.h>
#include <blaze/math/blas/trmv.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/MatVecMultExpr.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/expressions/VecScalarMultExpr.h>
#include <blaze/math/Functions.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsBlasCompatible.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/system/BLAS.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsComplexDouble.h>
#include <blaze/util/typetraits/IsComplexFloat.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATDVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// The TDMatDVecMultExpr class represents the compile time expression for multiplications
// between column-major dense matrices and dense vectors.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side dense vector
class TDMatDVecMultExpr : public DenseVector< TDMatDVecMultExpr<MT,VT>, false >
                        , private MatVecMultExpr
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     MRT;  //!< Result type of the left-hand side dense matrix expression.
   typedef typename VT::ResultType     VRT;  //!< Result type of the right-hand side dense vector expression.
   typedef typename MRT::ElementType   MET;  //!< Element type of the left-hand side dense matrix expression.
   typedef typename VRT::ElementType   VET;  //!< Element type of the right-hand side dense vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the left-hand side dense matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   enum { evaluateMatrix = ( IsComputation<MT>::value && IsSame<MET,VET>::value &&
                             IsBlasCompatible<MET>::value ) || RequiresEvaluation<MT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense vector expression.
   enum { evaluateVector = IsComputation<VT>::value || RequiresEvaluation<VT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The UseSMPAssign struct is a helper struct for the selection of the parallel evaluation
       strategy. In case either the matrix or the vector operand requires an intermediate
       evaluation, the nested \a value will be set to 1, otherwise it will be 0. */
   template< typename T1 >
   struct UseSMPAssign {
      enum { value = ( evaluateMatrix || evaluateVector ) };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the matrix type and the two involved vector types are suited for a BLAS kernel,
       the nested \a value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseBlasKernel {
      enum { value = BLAZE_BLAS_MODE &&
                     HasMutableDataAccess<T1>::value &&
                     HasConstDataAccess<T2>::value &&
                     HasConstDataAccess<T3>::value &&
                     !IsDiagonal<T2>::value &&
                     T1::vectorizable && T2::vectorizable && T3::vectorizable &&
                     IsBlasCompatible<typename T1::ElementType>::value &&
                     IsBlasCompatible<typename T2::ElementType>::value &&
                     IsBlasCompatible<typename T3::ElementType>::value &&
                     IsSame< typename T1::ElementType, typename T2::ElementType >::value &&
                     IsSame< typename T1::ElementType, typename T3::ElementType >::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the matrix type and the two involved vector types are suited for a vectorized
       computation of the matrix/vector multiplication, the nested \a value will be set to 1,
       otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseVectorizedDefaultKernel {
      enum { value = useOptimizedKernels &&
                     !IsDiagonal<T2>::value &&
                     T1::vectorizable && T2::vectorizable && T3::vectorizable &&
                     IsSame<typename T1::ElementType,typename T2::ElementType>::value &&
                     IsSame<typename T1::ElementType,typename T3::ElementType>::value &&
                     IntrinsicTrait<typename T1::ElementType>::addition &&
                     IntrinsicTrait<typename T1::ElementType>::multiplication };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TDMatDVecMultExpr<MT,VT>                    This;           //!< Type of this TDMatDVecMultExpr instance.
   typedef typename MultTrait<MRT,VRT>::Type           ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const ElementType                           ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                            CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side dense matrix operand.
   typedef typename SelectType< evaluateMatrix, const MRT, MCT >::Type  LT;

   //! Type for the assignment of the right-hand side dense vector operand.
   typedef typename SelectType< evaluateVector, const VRT, VCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = !IsDiagonal<MT>::value &&
                         MT::vectorizable && VT::vectorizable &&
                         IsSame<MET,VET>::value &&
                         IntrinsicTrait<MET>::addition &&
                         IntrinsicTrait<MET>::multiplication };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateMatrix && MT::smpAssignable &&
                          !evaluateVector && VT::smpAssignable };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDMatDVecMultExpr class.
   //
   // \param mat The left-hand side matrix operand of the multiplication expression.
   // \param vec The right-hand side vector operand of the multiplication expression.
   */
   explicit inline TDMatDVecMultExpr( const MT& mat, const VT& vec )
      : mat_( mat )  // Left-hand side dense matrix of the multiplication expression
      , vec_( vec )  // Right-hand side dense vector of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( mat_.columns() == vec_.size(), "Invalid matrix and vector sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < mat_.rows(), "Invalid vector access index" );

      if( ( IsStrictlyLower<MT>::value && index == 0UL ) ||
          ( IsStrictlyUpper<MT>::value && index == mat_.rows()-1UL ) ||
          mat_.columns() == 0UL )
         return ElementType();

      if( IsDiagonal<MT>::value )
         return mat_(index,index) * vec_[index];

      const size_t jbegin( ( IsUpper<MT>::value )
                           ?( IsStrictlyUpper<MT>::value ? index+1UL : index )
                           :( 0UL ) );
      const size_t jend( ( IsLower<MT>::value )
                         ?( IsStrictlyLower<MT>::value ? index : index+1UL )
                         :( mat_.columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jnum( jend - jbegin );
      const size_t jpos( jbegin + ( ( jnum - 1UL ) & size_t(-2) ) + 1UL );

      ElementType res( mat_(index,jbegin) * vec_[jbegin] );

      for( size_t j=jbegin+1UL; j<jpos; j+=2UL ) {
         res += mat_(index,j) * vec_[j] + mat_(index,j+1) * vec_[j+1UL];
      }
      if( jpos < jend ) {
         res += mat_(index,jpos) * vec_[jpos];
      }

      return res;
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
      if( index >= mat_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return mat_.rows();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side transpose dense matrix operand.
   //
   // \return The left-hand side transpose dense matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return mat_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const {
      return vec_;
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
      return ( mat_.isAliased( alias ) || vec_.isAliased( alias ) );
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
      return ( mat_.isAliased( alias ) || vec_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const {
      return mat_.isAligned() && vec_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return ( !BLAZE_BLAS_IS_PARALLEL ||
               ( IsComputation<MT>::value && !evaluateMatrix ) ||
               ( mat_.rows() * mat_.columns() < TDMATDVECMULT_THRESHOLD ) ) &&
             ( size() > SMP_TDMATDVECMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand vec_;  //!< Right-hand side dense vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-dense vector multiplication to a dense vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }
      else if( rhs.mat_.columns() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT A( serial( rhs.mat_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( rhs.vec_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      TDMatDVecMultExpr::selectAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense vectors (kernel selection)**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an assignment of a transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      if( ( IsDiagonal<MT1>::value ) ||
          ( IsComputation<MT>::value && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAssignKernel( y, A, x );
      else
         selectBlasAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default assignment kernel for the transpose dense matrix-dense
   // vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      if( IsStrictlyLower<MT1>::value ) {
         reset( y[0] );
      }

      if( !IsUpper<MT1>::value )
      {
         for( size_t i=( IsStrictlyLower<MT1>::value ? 1UL : 0UL ); i<M; ++i ) {
            y[i] = A(i,0UL) * x[0UL];
         }
      }

      for( size_t j=( IsUpper<MT1>::value && !IsStrictlyUpper<MT1>::value ? 0UL : 1UL ); j<N; ++j )
      {
         if( IsDiagonal<MT1>::value )
         {
            y[j] = A(j,j) * x[j];
         }
         else
         {
            const size_t ibegin( ( IsLower<MT1>::value )
                                 ?( IsStrictlyLower<MT1>::value ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( IsStrictlyUpper<MT1>::value ? j-1UL : j )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-2) ) );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] += A(i    ,j) * x[j];
               y[i+1UL] += A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] += A(ipos,j) * x[j];
            }
            if( IsUpper<MT1>::value ) {
               y[iend] = A(iend,j) * x[j];
            }
         }
      }

      if( IsStrictlyUpper<MT1>::value ) {
         reset( y[M-1UL] );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors (small matrices)****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a small transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (small matrices)*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the transpose dense
   // matrix-dense vector multiplication. This kernel is optimized for small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t ipos( remainder ? ( M & size_t(-IT::size) ) : M );
      BLAZE_INTERNAL_ASSERT( !remainder || ( M - ( M % IT::size ) ) == ipos, "Invalid end calculation" );

      size_t i( 0UL );

      for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*8UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
         y.store( i+IT::size*3UL, xmm4 );
         y.store( i+IT::size*4UL, xmm5 );
         y.store( i+IT::size*5UL, xmm6 );
         y.store( i+IT::size*6UL, xmm7 );
         y.store( i+IT::size*7UL, xmm8 );
      }

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*4UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
         y.store( i+IT::size*3UL, xmm4 );
      }

      for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*3UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
      }

      for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*2UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i         ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
         }

         y.store( i         , xmm1 );
         y.store( i+IT::size, xmm2 );
      }

      for( ; i<ipos; i+=IT::size )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1;

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 = xmm1 + A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1 );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+1UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         ElementType value = ElementType();

         for( size_t j=jbegin; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] = value;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors (large matrices)****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a large transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (large matrices)*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the transpose dense
   // matrix-dense vector multiplication. This kernel is optimized for large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % IT::size ) == 0UL, "Invalid block size detected" );

      reset( y );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( min( itmp, ( IsStrictlyUpper<MT1>::value ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? ( iend & size_t(-IT::size) ) : iend );
            BLAZE_INTERNAL_ASSERT( !remainder || ( iend - ( iend % IT::size ) ) == ipos, "Invalid end calculation" );

            size_t i( ( IsLower<MT1>::value )
                      ?( max( ii, ( IsStrictlyLower<MT1>::value ? jj+1UL : jj ) & size_t(-IT::size) ) )
                      :( ii ) );

            for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
                  xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
                  xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
                  xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
                  xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3 );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4 );
               y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) + xmm5 );
               y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) + xmm6 );
               y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) + xmm7 );
               y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) + xmm8 );
            }

            for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3 );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4 );
            }

            for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
            {
               IntrinsicType xmm1, xmm2, xmm3;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3 );
            }

            for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
            {
               IntrinsicType xmm1, xmm2;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i         ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1 );
               y.store( i+IT::size, y.load(i+IT::size) + xmm2 );
            }

            for( ; i<ipos; i+=IT::size )
            {
               IntrinsicType xmm1;

               for( size_t j=jj; j<jend; ++j ) {
                  xmm1 = xmm1 + A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1 );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value = ElementType();

               for( size_t j=jj; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (default)********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a large transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseBlasKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectLargeAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors******************************************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication based on the
   // according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseBlasKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef typename VT1::ElementType  ET;

      if( IsTriangular<MT1>::value ) {
         assign( y, x );
         trmv( y, A, ( IsLower<MT1>::value )?( CblasLower ):( CblasUpper ) );
      }
      else {
         gemv( y, A, x, ET(1), ET(0) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-dense vector multiplication to a sparse vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ) {
         return;
      }

      LT A( serial( rhs.mat_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( rhs.vec_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      TDMatDVecMultExpr::selectAddAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors (kernel selection)*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an addition assignment of a transpose dense matrix-dense
   //        vector multiplication to a dense vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      if( ( IsDiagonal<MT1>::value ) ||
          ( IsComputation<MT>::value && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAddAssignKernel( y, A, x );
      else
         selectBlasAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default addition assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         if( IsDiagonal<MT1>::value )
         {
            y[j] += A(j,j) * x[j];
         }
         else
         {
            const size_t ibegin( ( IsLower<MT1>::value )
                                 ?( IsStrictlyLower<MT1>::value ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( IsStrictlyUpper<MT1>::value ? j : j+1UL )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-2) ) );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] += A(i    ,j) * x[j];
               y[i+1UL] += A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] += A(ipos,j) * x[j];
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (small matrices)*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (small matrices)********************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t ipos( remainder ? ( M & size_t(-IT::size) ) : M );
      BLAZE_INTERNAL_ASSERT( !remainder || ( M - ( M % IT::size ) ) == ipos, "Invalid end calculation" );

      size_t i( 0UL );

      for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*8UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i             ) );
         IntrinsicType xmm2( y.load(i+IT::size    ) );
         IntrinsicType xmm3( y.load(i+IT::size*2UL) );
         IntrinsicType xmm4( y.load(i+IT::size*3UL) );
         IntrinsicType xmm5( y.load(i+IT::size*4UL) );
         IntrinsicType xmm6( y.load(i+IT::size*5UL) );
         IntrinsicType xmm7( y.load(i+IT::size*6UL) );
         IntrinsicType xmm8( y.load(i+IT::size*7UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
         y.store( i+IT::size*3UL, xmm4 );
         y.store( i+IT::size*4UL, xmm5 );
         y.store( i+IT::size*5UL, xmm6 );
         y.store( i+IT::size*6UL, xmm7 );
         y.store( i+IT::size*7UL, xmm8 );
      }

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*4UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i             ) );
         IntrinsicType xmm2( y.load(i+IT::size    ) );
         IntrinsicType xmm3( y.load(i+IT::size*2UL) );
         IntrinsicType xmm4( y.load(i+IT::size*3UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
         y.store( i+IT::size*3UL, xmm4 );
      }

      for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*3UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i             ) );
         IntrinsicType xmm2( y.load(i+IT::size    ) );
         IntrinsicType xmm3( y.load(i+IT::size*2UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
      }

      for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*2UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i         ) );
         IntrinsicType xmm2( y.load(i+IT::size) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i         ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
         }

         y.store( i         , xmm1 );
         y.store( i+IT::size, xmm2 );
      }

      for( ; i<ipos; i+=IT::size )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i) );

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 = xmm1 + A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1 );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+1UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         ElementType value = ElementType();

         for( size_t j=jbegin; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] += value;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (large matrices)*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (large matrices)********************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % IT::size ) == 0UL, "Invalid block size detected" );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( min( itmp, ( IsStrictlyUpper<MT1>::value ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? ( iend & size_t(-IT::size) ) : iend );
            BLAZE_INTERNAL_ASSERT( !remainder || ( iend - ( iend % IT::size ) ) == ipos, "Invalid end calculation" );

            size_t i( ( IsLower<MT1>::value )
                      ?( max( ii, ( IsStrictlyLower<MT1>::value ? jj+1UL : jj ) & size_t(-IT::size) ) )
                      :( ii ) );

            for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
                  xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
                  xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
                  xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
                  xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3 );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4 );
               y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) + xmm5 );
               y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) + xmm6 );
               y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) + xmm7 );
               y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) + xmm8 );
            }

            for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3 );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4 );
            }

            for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
            {
               IntrinsicType xmm1, xmm2, xmm3;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3 );
            }

            for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
            {
               IntrinsicType xmm1, xmm2;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i         ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1 );
               y.store( i+IT::size, y.load(i+IT::size) + xmm2 );
            }

            for( ; i<ipos; i+=IT::size )
            {
               IntrinsicType xmm1;

               for( size_t j=jj; j<jend; ++j ) {
                  xmm1 = xmm1 + A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1 );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value = ElementType();

               for( size_t j=jj; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (default)***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a large
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseBlasKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectLargeAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors*********************************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose matrix-vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication based on the
   // according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseBlasKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef typename VT1::ElementType  ET;

      if( IsTriangular<MT1>::value ) {
         typename VT1::ResultType tmp( serial( x ) );
         trmv( tmp, A, ( IsLower<MT1>::value )?( CblasLower ):( CblasUpper ) );
         addAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(1), ET(1) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ) {
         return;
      }

      LT A( serial( rhs.mat_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( rhs.vec_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      TDMatDVecMultExpr::selectSubAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors (kernel selection)**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for a subtraction assignment of a transpose dense matrix-
   //        dense vector multiplication to a dense vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      if( ( IsDiagonal<MT1>::value ) ||
          ( IsComputation<MT>::value && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallSubAssignKernel( y, A, x );
      else
         selectBlasSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         if( IsDiagonal<MT1>::value )
         {
            y[j] -= A(j,j) * x[j];
         }
         else
         {
            const size_t ibegin( ( IsLower<MT1>::value )
                                 ?( IsStrictlyLower<MT1>::value ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( IsStrictlyUpper<MT1>::value ? j : j+1UL )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-2) ) );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] -= A(i    ,j) * x[j];
               y[i+1UL] -= A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] -= A(ipos,j) * x[j];
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (small matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (small matrices)*****************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a small transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for small
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t ipos( remainder ? ( M & size_t(-IT::size) ) : M );
      BLAZE_INTERNAL_ASSERT( !remainder || ( M - ( M % IT::size ) ) == ipos, "Invalid end calculation" );

      size_t i( 0UL );

      for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*8UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i             ) );
         IntrinsicType xmm2( y.load(i+IT::size    ) );
         IntrinsicType xmm3( y.load(i+IT::size*2UL) );
         IntrinsicType xmm4( y.load(i+IT::size*3UL) );
         IntrinsicType xmm5( y.load(i+IT::size*4UL) );
         IntrinsicType xmm6( y.load(i+IT::size*5UL) );
         IntrinsicType xmm7( y.load(i+IT::size*6UL) );
         IntrinsicType xmm8( y.load(i+IT::size*7UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.load(i             ,j) * x1;
            xmm2 = xmm2 - A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 - A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 - A.load(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 - A.load(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 - A.load(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 - A.load(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 - A.load(i+IT::size*7UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
         y.store( i+IT::size*3UL, xmm4 );
         y.store( i+IT::size*4UL, xmm5 );
         y.store( i+IT::size*5UL, xmm6 );
         y.store( i+IT::size*6UL, xmm7 );
         y.store( i+IT::size*7UL, xmm8 );
      }

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*4UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i             ) );
         IntrinsicType xmm2( y.load(i+IT::size    ) );
         IntrinsicType xmm3( y.load(i+IT::size*2UL) );
         IntrinsicType xmm4( y.load(i+IT::size*3UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.load(i             ,j) * x1;
            xmm2 = xmm2 - A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 - A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 - A.load(i+IT::size*3UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
         y.store( i+IT::size*3UL, xmm4 );
      }

      for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*3UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i             ) );
         IntrinsicType xmm2( y.load(i+IT::size    ) );
         IntrinsicType xmm3( y.load(i+IT::size*2UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.load(i             ,j) * x1;
            xmm2 = xmm2 - A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 - A.load(i+IT::size*2UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+IT::size    , xmm2 );
         y.store( i+IT::size*2UL, xmm3 );
      }

      for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*2UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i         ) );
         IntrinsicType xmm2( y.load(i+IT::size) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.load(i         ,j) * x1;
            xmm2 = xmm2 - A.load(i+IT::size,j) * x1;
         }

         y.store( i         , xmm1 );
         y.store( i+IT::size, xmm2 );
      }

      for( ; i<ipos; i+=IT::size )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1( y.load(i) );

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 = xmm1 - A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1 );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+1UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         ElementType value = ElementType();

         for( size_t j=jbegin; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] -= value;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (large matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (large matrices)*****************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a large transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for large
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % IT::size ) == 0UL, "Invalid block size detected" );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( min( itmp, ( IsStrictlyUpper<MT1>::value ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? ( iend & size_t(-IT::size) ) : iend );
            BLAZE_INTERNAL_ASSERT( !remainder || ( iend - ( iend % IT::size ) ) == ipos, "Invalid end calculation" );

            size_t i( ( IsLower<MT1>::value )
                      ?( max( ii, ( IsStrictlyLower<MT1>::value ? jj+1UL : jj ) & size_t(-IT::size) ) )
                      :( ii ) );

            for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
                  xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
                  xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
                  xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
                  xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3 );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) - xmm4 );
               y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) - xmm5 );
               y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) - xmm6 );
               y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) - xmm7 );
               y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) - xmm8 );
            }

            for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3 );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) - xmm4 );
            }

            for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
            {
               IntrinsicType xmm1, xmm2, xmm3;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1 );
               y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2 );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3 );
            }

            for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
            {
               IntrinsicType xmm1, xmm2;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i         ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
               }

               y.store( i         , y.load(i         ) - xmm1 );
               y.store( i+IT::size, y.load(i+IT::size) - xmm2 );
            }

            for( ; i<ipos; i+=IT::size )
            {
               IntrinsicType xmm1;

               for( size_t j=jj; j<jend; ++j ) {
                  xmm1 = xmm1 + A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) - xmm1 );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value = ElementType();

               for( size_t j=jj; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] -= value;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (default)********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a large
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseBlasKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectLargeSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors******************************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subtraction assignment of a transpose matrix-vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication based on the
   // according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseBlasKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef typename VT1::ElementType  ET;

      if( IsTriangular<MT1>::value ) {
         typename VT1::ResultType tmp( serial( x ) );
         trmv( tmp, A, ( IsLower<MT1>::value )?( CblasLower ):( CblasUpper ) );
         subAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(-1), ET(1) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a transpose dense matrix-dense vector multiplication to
   //        a dense vector (\f$ \vec{y}*=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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
   /*!\brief SMP assignment of a transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }
      else if( rhs.mat_.columns() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      smpAssign( ~lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-dense vector multiplication to a sparse
   //        vector (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-dense vector multiplication expression to a sparse vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpAssign( SparseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpAddAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ) {
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      smpAddAssign( ~lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpSubAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ) {
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      smpSubAssign( ~lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}*=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // transpose dense matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpMultAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATVECMULTEXPR( MT, VT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DVECSCALARMULTEXPR SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expression object for scaled transpose dense matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// This specialization of the DVecScalarMultExpr class represents the compile time expression
// for scaled multiplications between a column-major dense matrix and a non-transpose dense
// vector.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT    // Type of the right-hand side dense vector
        , typename ST >  // Type of the side scalar value
class DVecScalarMultExpr< TDMatDVecMultExpr<MT,VT>, ST, false >
   : public DenseVector< DVecScalarMultExpr< TDMatDVecMultExpr<MT,VT>, ST, false >, false >
   , private VecScalarMultExpr
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef TDMatDVecMultExpr<MT,VT>    MVM;  //!< Type of the transpose dense matrix-dense vector multiplication expression.
   typedef typename MVM::ResultType    RES;  //!< Result type of the dense matrix-dense vector multiplication expression.
   typedef typename MT::ResultType     MRT;  //!< Result type of the left-hand side dense matrix expression.
   typedef typename VT::ResultType     VRT;  //!< Result type of the right-hand side dense vector expression.
   typedef typename MRT::ElementType   MET;  //!< Element type of the left-hand side dense matrix expression.
   typedef typename VRT::ElementType   VET;  //!< Element type of the right-hand side dense vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the left-hand side dense matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   enum { evaluateMatrix = ( IsComputation<MT>::value && IsSame<MET,VET>::value &&
                             IsBlasCompatible<MET>::value ) || RequiresEvaluation<MT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense vector expression.
   enum { evaluateVector = IsComputation<VT>::value || RequiresEvaluation<VT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The UseSMPAssign struct is a helper struct for the selection of the parallel evaluation
       strategy. In case either the matrix or the vector operand requires an intermediate
       evaluation, the nested \a value will be set to 1, otherwise it will be 0. */
   template< typename T1 >
   struct UseSMPAssign {
      enum { value = ( evaluateMatrix || evaluateVector ) };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the matrix type, the two involved vector types, and the scalar type are suited
       for a BLAS kernel, the nested \a value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   struct UseBlasKernel {
      enum { value = BLAZE_BLAS_MODE &&
                     HasMutableDataAccess<T1>::value &&
                     HasConstDataAccess<T2>::value &&
                     HasConstDataAccess<T3>::value &&
                     !IsDiagonal<T2>::value &&
                     T1::vectorizable && T2::vectorizable && T3::vectorizable &&
                     IsBlasCompatible<typename T1::ElementType>::value &&
                     IsBlasCompatible<typename T2::ElementType>::value &&
                     IsBlasCompatible<typename T3::ElementType>::value &&
                     IsSame< typename T1::ElementType, typename T2::ElementType >::value &&
                     IsSame< typename T1::ElementType, typename T3::ElementType >::value &&
                     !( IsBuiltin<typename T1::ElementType>::value && IsComplex<T4>::value ) };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the two involved vector types, the matrix type, and the scalar type are suited
       for a vectorized computation of the scaled vector/matrix multiplication, the nested
       \a value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   struct UseVectorizedDefaultKernel {
      enum { value = useOptimizedKernels &&
                     !IsDiagonal<T2>::value &&
                     T1::vectorizable && T2::vectorizable && T3::vectorizable &&
                     IsSame<typename T1::ElementType,typename T2::ElementType>::value &&
                     IsSame<typename T1::ElementType,typename T3::ElementType>::value &&
                     IsSame<typename T1::ElementType,T4>::value &&
                     IntrinsicTrait<typename T1::ElementType>::addition &&
                     IntrinsicTrait<typename T1::ElementType>::multiplication };
   };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DVecScalarMultExpr<MVM,ST,false>            This;           //!< Type of this DVecScalarMultExpr instance.
   typedef typename MultTrait<RES,ST>::Type            ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const ElementType                           ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                            CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense vector expression.
   typedef const TDMatDVecMultExpr<MT,VT>  LeftOperand;

   //! Composite type of the right-hand side scalar value.
   typedef ST  RightOperand;

   //! Type for the assignment of the dense matrix operand of the left-hand side expression.
   typedef typename SelectType< evaluateMatrix, const MRT, MCT >::Type  LT;

   //! Type for the assignment of the dense vector operand of the left-hand side expression.
   typedef typename SelectType< evaluateVector, const VRT, VCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = !IsDiagonal<MT>::value &&
                         MT::vectorizable && VT::vectorizable &&
                         IsSame<MET,VET>::value &&
                         IsSame<MET,ST>::value &&
                         IntrinsicTrait<MET>::addition &&
                         IntrinsicTrait<MET>::multiplication };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateMatrix && MT::smpAssignable &&
                          !evaluateVector && VT::smpAssignable };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecScalarMultExpr class.
   //
   // \param vector The left-hand side dense vector of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   explicit inline DVecScalarMultExpr( const MVM& vector, ST scalar )
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
      return vector_.canAlias( alias );
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
      typename MVM::LeftOperand A( vector_.leftOperand() );
      return ( !BLAZE_BLAS_IS_PARALLEL ||
               ( IsComputation<MT>::value && !evaluateMatrix ) ||
               ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) ) &&
             ( size() > SMP_TDMATDVECMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vector_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( right ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      DVecScalarMultExpr::selectAssignKernel( ~lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Assignment to dense vectors (kernel selection)**********************************************
   /*!\brief Selection of the kernel for an assignment of a scaled transpose dense matrix-dense
   //        vector multiplication to a dense vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      if( ( IsDiagonal<MT1>::value ) ||
          ( IsComputation<MT>::value && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAssignKernel( y, A, x, scalar );
      else
         selectBlasAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the default assignment kernel for the scaled transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      if( IsStrictlyLower<MT1>::value ) {
         reset( y[0] );
      }

      if( !IsUpper<MT1>::value )
      {
         for( size_t i=( IsStrictlyLower<MT1>::value ? 1UL : 0UL ); i<M; ++i ) {
            y[i] = A(i,0UL) * x[0UL];
         }
      }

      for( size_t j=( IsUpper<MT1>::value && !IsStrictlyUpper<MT1>::value ? 0UL : 1UL ); j<N; ++j )
      {
         if( IsDiagonal<MT1>::value )
         {
            y[j] = A(j,j) * x[j] * scalar;
         }
         else
         {
            const size_t ibegin( ( IsLower<MT1>::value )
                                 ?( IsStrictlyLower<MT1>::value ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( IsStrictlyUpper<MT1>::value ? j-1UL : j )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-2) ) );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] += A(i    ,j) * x[j];
               y[i+1UL] += A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] += A(ipos,j) * x[j];
            }
            if( IsUpper<MT1>::value ) {
               y[iend] = A(iend,j) * x[j];
            }
         }
      }

      if( IsStrictlyUpper<MT1>::value ) {
         reset( y[M-1UL] );
      }

      if( !IsDiagonal<MT1>::value )
      {
         const size_t iend( IsStrictlyUpper<MT1>::value ? M-1UL : M );
         for( size_t i=( IsStrictlyLower<MT1>::value ? 1UL : 0UL ); i<iend; ++i ) {
            y[i] *= scalar;
         }
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors (small matrices)****************************************
   /*!\brief Default assignment of a small scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (small matrices)*****************************
   /*!\brief Vectorized default assignment of a small scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t ipos( remainder ? ( M & size_t(-IT::size) ) : M );
      BLAZE_INTERNAL_ASSERT( !remainder || ( M - ( M % IT::size ) ) == ipos, "Invalid end calculation" );

      const IntrinsicType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*8UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
         }

         y.store( i             , xmm1*factor );
         y.store( i+IT::size    , xmm2*factor );
         y.store( i+IT::size*2UL, xmm3*factor );
         y.store( i+IT::size*3UL, xmm4*factor );
         y.store( i+IT::size*4UL, xmm5*factor );
         y.store( i+IT::size*5UL, xmm6*factor );
         y.store( i+IT::size*6UL, xmm7*factor );
         y.store( i+IT::size*7UL, xmm8*factor );
      }

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*4UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
         }

         y.store( i             , xmm1*factor );
         y.store( i+IT::size    , xmm2*factor );
         y.store( i+IT::size*2UL, xmm3*factor );
         y.store( i+IT::size*3UL, xmm4*factor );
      }

      for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*3UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
         }

         y.store( i             , xmm1*factor );
         y.store( i+IT::size    , xmm2*factor );
         y.store( i+IT::size*2UL, xmm3*factor );
      }

      for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*2UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i         ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
         }

         y.store( i         , xmm1*factor );
         y.store( i+IT::size, xmm2*factor );
      }

      for( ; i<ipos; i+=IT::size )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i,j) * x1;
         }

         y.store( i, xmm1*factor );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+1UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         ElementType value = ElementType();

         for( size_t j=jbegin; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] = value * scalar;
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors (large matrices)****************************************
   /*!\brief Default assignment of a large scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (large matrices)*****************************
   /*!\brief Vectorized default assignment of a large scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % IT::size ) == 0UL, "Invalid block size detected" );

      const IntrinsicType factor( set( scalar ) );

      reset( y );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( min( itmp, ( IsStrictlyUpper<MT1>::value ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? ( iend & size_t(-IT::size) ) : iend );
            BLAZE_INTERNAL_ASSERT( !remainder || ( iend - ( iend % IT::size ) ) == ipos, "Invalid end calculation" );

            size_t i( ( IsLower<MT1>::value )
                      ?( max( ii, ( IsStrictlyLower<MT1>::value ? jj+1UL : jj ) & size_t(-IT::size) ) )
                      :( ii ) );

            for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
                  xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
                  xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
                  xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
                  xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4*factor );
               y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) + xmm5*factor );
               y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) + xmm6*factor );
               y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) + xmm7*factor );
               y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) + xmm8*factor );
            }

            for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4*factor );
            }

            for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
            {
               IntrinsicType xmm1, xmm2, xmm3;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
            }

            for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
            {
               IntrinsicType xmm1, xmm2;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i         ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1*factor );
               y.store( i+IT::size, y.load(i+IT::size) + xmm2*factor );
            }

            for( ; i<ipos; i+=IT::size )
            {
               IntrinsicType xmm1;

               for( size_t j=jj; j<jend; ++j ) {
                  xmm1 = xmm1 + A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1*factor );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value = ElementType();

               for( size_t j=jj; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (default)********************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a large scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseBlasKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectLargeAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors******************************************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication based
   // on the according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseBlasKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef typename VT1::ElementType  ET;

      if( IsTriangular<MT1>::value ) {
         assign( y, scalar * x );
         trmv( y, A, ( IsLower<MT1>::value )?( CblasLower ):( CblasUpper ) );
      }
      else {
         gemv( y, A, x, ET(scalar), ET(0) );
      }
   }
#endif
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-dense vector multiplication to a sparse
   //        vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( right ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      DVecScalarMultExpr::selectAddAssignKernel( ~lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors (kernel selection)*************************************
   /*!\brief Selection of the kernel for an addition assignment of a scaled transpose dense
   //        matrix-dense vector multiplication to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      if( ( IsDiagonal<MT1>::value ) ||
          ( IsComputation<MT>::value && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAddAssignKernel( y, A, x, scalar );
      else
         selectBlasAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      y.addAssign( A * x * scalar );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (small matrices)*******************************
   /*!\brief Default addition assignment of a small scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (small matrices)********************
   /*!\brief Vectorized default addition assignment of a small scaled transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the scaled
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for small
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t ipos( remainder ? ( M & size_t(-IT::size) ) : M );
      BLAZE_INTERNAL_ASSERT( !remainder || ( M - ( M % IT::size ) ) == ipos, "Invalid end calculation" );

      const IntrinsicType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*8UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
         }

         y.store( i             , y.load(i             ) + xmm1*factor );
         y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
         y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
         y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4*factor );
         y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) + xmm5*factor );
         y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) + xmm6*factor );
         y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) + xmm7*factor );
         y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) + xmm8*factor );
      }

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*4UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
         }

         y.store( i             , y.load(i             ) + xmm1*factor );
         y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
         y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
         y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4*factor );
      }

      for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*3UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
         }

         y.store( i             , y.load(i             ) + xmm1*factor );
         y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
         y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
      }

      for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*2UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i         ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
         }

         y.store( i         , y.load(i         ) + xmm1*factor );
         y.store( i+IT::size, y.load(i+IT::size) + xmm2*factor );
      }

      for( ; i<ipos; i+=IT::size )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1;

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 = xmm1 + A.load(i,j) * set( x[j] );
         }

         y.store( i, y.load(i) + xmm1*factor );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+1UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         ElementType value = ElementType();

         for( size_t j=jbegin; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] += value * scalar;
      }
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (large matrices)*******************************
   /*!\brief Default addition assignment of a large scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (large matrices)********************
   /*!\brief Vectorized default addition assignment of a large scaled transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the scaled
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for large
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % IT::size ) == 0UL, "Invalid block size detected" );

      const IntrinsicType factor( set( scalar ) );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( min( itmp, ( IsStrictlyUpper<MT1>::value ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? ( iend & size_t(-IT::size) ) : iend );
            BLAZE_INTERNAL_ASSERT( !remainder || ( iend - ( iend % IT::size ) ) == ipos, "Invalid end calculation" );

            size_t i( ( IsLower<MT1>::value )
                      ?( max( ii, ( IsStrictlyLower<MT1>::value ? jj+1UL : jj ) & size_t(-IT::size) ) )
                      :( ii ) );

            for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
                  xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
                  xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
                  xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
                  xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4*factor );
               y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) + xmm5*factor );
               y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) + xmm6*factor );
               y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) + xmm7*factor );
               y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) + xmm8*factor );
            }

            for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) + xmm4*factor );
            }

            for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
            {
               IntrinsicType xmm1, xmm2, xmm3;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) + xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) + xmm3*factor );
            }

            for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
            {
               IntrinsicType xmm1, xmm2;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i         ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1*factor );
               y.store( i+IT::size, y.load(i+IT::size) + xmm2*factor );
            }

            for( ; i<ipos; i+=IT::size )
            {
               IntrinsicType xmm1;

               for( size_t j=jj; j<jend; ++j ) {
                  xmm1 = xmm1 + A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1*factor );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value = ElementType();

               for( size_t j=jj; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (default)***********************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a large
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseBlasKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectLargeAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors*********************************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication based
   // on the according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseBlasKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef typename VT1::ElementType  ET;

      if( IsTriangular<MT1>::value ) {
         typename VT1::ResultType tmp( serial( scalar * x ) );
         trmv( tmp, A, ( IsLower<MT1>::value )?( CblasLower ):( CblasUpper ) );
         addAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(scalar), ET(1) );
      }
   }
#endif
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a scaled
   // dense transpose matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( right ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      DVecScalarMultExpr::selectSubAssignKernel( ~lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors (kernel selection)**********************************
   /*!\brief Selection of the kernel for a subtraction assignment of a scaled transpose dense
   //        matrix-dense vector multiplication to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      if( ( IsDiagonal<MT1>::value ) ||
          ( IsComputation<MT>::value && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallSubAssignKernel( y, A, x, scalar );
      else
         selectBlasSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      y.subAssign( A * x * scalar );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (small matrices)****************************
   /*!\brief Default subtraction assignment of a small scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (small matrices)*****************
   /*!\brief Vectorized default subtraction assignment of a small scaled transpose dense matrix-
   //        dense vector multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // scaled transpose dense matrix-dense vector multiplication. This kernel is optimized for
   // small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t ipos( remainder ? ( M & size_t(-IT::size) ) : M );
      BLAZE_INTERNAL_ASSERT( !remainder || ( M - ( M % IT::size ) ) == ipos, "Invalid end calculation" );

      const IntrinsicType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*8UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
         }

         y.store( i             , y.load(i             ) - xmm1*factor );
         y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2*factor );
         y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3*factor );
         y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) - xmm4*factor );
         y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) - xmm5*factor );
         y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) - xmm6*factor );
         y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) - xmm7*factor );
         y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) - xmm8*factor );
      }

      for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*4UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3, xmm4;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
         }

         y.store( i             , y.load(i             ) - xmm1*factor );
         y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2*factor );
         y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3*factor );
         y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) - xmm4*factor );
      }

      for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*3UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2, xmm3;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i             ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
         }

         y.store( i             , y.load(i             ) - xmm1*factor );
         y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2*factor );
         y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3*factor );
      }

      for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size*2UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1, xmm2;

         for( size_t j=jbegin; j<jend; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.load(i         ,j) * x1;
            xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
         }

         y.store( i         , y.load(i         ) - xmm1*factor );
         y.store( i+IT::size, y.load(i+IT::size) - xmm2*factor );
      }

      for( ; i<ipos; i+=IT::size )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+IT::size, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         IntrinsicType xmm1;

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 = xmm1 + A.load(i,j) * set( x[j] );
         }

         y.store( i, y.load(i) - xmm1*factor );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper<MT1>::value )
                              ?( IsStrictlyUpper<MT1>::value ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower<MT1>::value )
                            ?( min( i+1UL, N ) - ( IsStrictlyLower<MT1>::value ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         ElementType value = ElementType();

         for( size_t j=jbegin; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] -= value * scalar;
      }
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (large matrices)****************************
   /*!\brief Default subtraction assignment of a large scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (large matrices)*****************
   /*!\brief Vectorized default subtraction assignment of a large scaled transpose dense matrix-
   //        dense vector multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // scaled transpose dense matrix-dense vector multiplication. This kernel is optimized for
   // large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const bool remainder( !IsPadded<MT1>::value || !IsPadded<VT1>::value );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % IT::size ) == 0UL, "Invalid block size detected" );

      const IntrinsicType factor( set( scalar ) );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper<MT1>::value )
                               ?( min( itmp, ( IsStrictlyUpper<MT1>::value ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? ( iend & size_t(-IT::size) ) : iend );
            BLAZE_INTERNAL_ASSERT( !remainder || ( iend - ( iend % IT::size ) ) == ipos, "Invalid end calculation" );

            size_t i( ( IsLower<MT1>::value )
                      ?( max( ii, ( IsStrictlyLower<MT1>::value ? jj+1UL : jj ) & size_t(-IT::size) ) )
                      :( ii ) );

            for( ; (i+IT::size*7UL) < ipos; i+=IT::size*8UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
                  xmm5 = xmm5 + A.load(i+IT::size*4UL,j) * x1;
                  xmm6 = xmm6 + A.load(i+IT::size*5UL,j) * x1;
                  xmm7 = xmm7 + A.load(i+IT::size*6UL,j) * x1;
                  xmm8 = xmm8 + A.load(i+IT::size*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3*factor );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) - xmm4*factor );
               y.store( i+IT::size*4UL, y.load(i+IT::size*4UL) - xmm5*factor );
               y.store( i+IT::size*5UL, y.load(i+IT::size*5UL) - xmm6*factor );
               y.store( i+IT::size*6UL, y.load(i+IT::size*6UL) - xmm7*factor );
               y.store( i+IT::size*7UL, y.load(i+IT::size*7UL) - xmm8*factor );
            }

            for( ; (i+IT::size*3UL) < ipos; i+=IT::size*4UL )
            {
               IntrinsicType xmm1, xmm2, xmm3, xmm4;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
                  xmm4 = xmm4 + A.load(i+IT::size*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3*factor );
               y.store( i+IT::size*3UL, y.load(i+IT::size*3UL) - xmm4*factor );
            }

            for( ; (i+IT::size*2UL) < ipos; i+=IT::size*3UL )
            {
               IntrinsicType xmm1, xmm2, xmm3;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i             ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size    ,j) * x1;
                  xmm3 = xmm3 + A.load(i+IT::size*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1*factor );
               y.store( i+IT::size    , y.load(i+IT::size    ) - xmm2*factor );
               y.store( i+IT::size*2UL, y.load(i+IT::size*2UL) - xmm3*factor );
            }

            for( ; (i+IT::size) < ipos; i+=IT::size*2UL )
            {
               IntrinsicType xmm1, xmm2;

               for( size_t j=jj; j<jend; ++j ) {
                  const IntrinsicType x1( set( x[j] ) );
                  xmm1 = xmm1 + A.load(i         ,j) * x1;
                  xmm2 = xmm2 + A.load(i+IT::size,j) * x1;
               }

               y.store( i         , y.load(i         ) - xmm1*factor );
               y.store( i+IT::size, y.load(i+IT::size) - xmm2*factor );
            }

            for( ; i<ipos; i+=IT::size )
            {
               IntrinsicType xmm1;

               for( size_t j=jj; j<jend; ++j ) {
                  xmm1 = xmm1 + A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) - xmm1*factor );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value = ElementType();

               for( size_t j=jj; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] -= value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (default)********************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a large
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseBlasKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectLargeSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors******************************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication based
   // on the according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseBlasKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef typename VT1::ElementType  ET;

      if( IsTriangular<MT1>::value ) {
         typename VT1::ResultType tmp( serial( scalar * x ) );
         trmv( tmp, A, ( IsLower<MT1>::value )?( CblasLower ):( CblasUpper ) );
         subAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(-scalar), ET(1) );
      }
   }
#endif
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}*=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors**************************************************************
   /*!\brief SMP assignment of a scaled transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      smpAssign( ~lhs, A * x * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*!\brief SMP assignment of a scaled transpose dense matrix-dense vector multiplication to a
   //        sparse vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a sparse vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpAssign( SparseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*!\brief SMP addition assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpAddAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      smpAddAssign( ~lhs, A * x * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*!\brief SMP subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a scaled
   // dense transpose matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpSubAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      smpSubAssign( ~lhs, A * x * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*!\brief SMP multiplication assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}*=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   // Due to the explicit application of the SFINAE principle, this function can only be
   // selected by the compiler in case the expression specific parallel evaluation strategy
   // is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT1> >::Type
      smpMultAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( MVM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( MVM );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a column-major dense matrix and a dense
//        vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side column-major dense matrix for the multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a column-major dense matrix and a dense
// vector:

   \code
   using blaze::columnMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<double,columnMajor> A;
   blaze::DynamicVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved element types \a T1::ElementType and \a T2::ElementType. Both the
// dense matrix type \a T1 and the dense vector type \a T2 as well as the two element types
// \a T1::ElementType and \a T2::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side dense vector
inline const typename DisableIf< IsMatMatMultExpr<T1>, TDMatDVecMultExpr<T1,T2> >::Type
   operator*( const DenseMatrix<T1,true>& mat, const DenseVector<T2,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   if( (~mat).columns() != (~vec).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix and vector sizes do not match" );
   }

   return TDMatDVecMultExpr<T1,T2>( ~mat, ~vec );
}
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT >
struct Size< TDMatDVecMultExpr<MT,VT> > : public Rows<MT>
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
template< typename MT, typename VT >
struct IsAligned< TDMatDVecMultExpr<MT,VT> >
   : public IsTrue< And< IsAligned<MT>, IsAligned<VT> >::value >
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
template< typename MT, typename VT, bool AF >
struct SubvectorExprTrait< TDMatDVecMultExpr<MT,VT>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename MultExprTrait< typename SubmatrixExprTrait<const MT,AF>::Type
                                 , typename SubvectorExprTrait<const VT,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
