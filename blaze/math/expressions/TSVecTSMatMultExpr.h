//=================================================================================================
/*!
//  \file blaze/math/expressions/TSVecTSMatMultExpr.h
//  \brief Header file for the transpose sparse vector/transpose sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TSVECTSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSVECTSMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TSVECTSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse vector-sparse matrix multiplications.
// \ingroup sparse_vector_expression
//
// The TSVecTSMatMultExpr class represents the compile time expression for multiplications
// between transpose sparse vectors and column-major sparse matrices.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
class TSVecTSMatMultExpr : public SparseVector< TSVecTSMatMultExpr<VT,MT>, true >
                         , private TVecMatMultExpr
                         , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     VRT;  //!< Result type of the left-hand side sparse vector expression.
   typedef typename MT::ResultType     MRT;  //!< Result type of the right-hand side sparse matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the left-hand side sparse vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side sparse vector expression.
   enum { evaluateVector = IsComputation<VT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse matrix expression.
   enum { evaluateMatrix = RequiresEvaluation<MT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the either the vector or the matrix operand require an intermediate evaluation,
       the nested \value will be set to 0, otherwise it will be 1. */
   template< typename T1, typename T2, typename T3 >
   struct UseSMPAssignKernel {
      enum { value = evaluateVector || evaluateMatrix };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TSVecTSMatMultExpr<VT,MT>           This;           //!< Type of this TSVecTSMatMultExpr instance.
   typedef typename MultTrait<VRT,MRT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.
   typedef const ElementType                   ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                    CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side sparse vector operand.
   typedef typename SelectType< evaluateVector, const VRT, VCT >::Type  LT;

   //! Type for the assignment of the right-hand side sparse matrix operand.
   typedef MCT  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateMatrix && !evaluateVector };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TSVecTSMatMultExpr class.
   //
   // \param vec The left-hand side sparse vector operand of the multiplication expression.
   // \param mat The right-hand side sparse matrix operand of the multiplication expression.
   */
   explicit inline TSVecTSMatMultExpr( const VT& vec, const MT& mat )
      : vec_( vec )  // Left-hand side sparse vector of the multiplication expression
      , mat_( mat )  // Right-hand side sparse matrix of the multiplication expression
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

      typedef typename RemoveReference<VCT>::Type::ConstIterator  VectorIterator;
      typedef typename RemoveReference<MCT>::Type::ConstIterator  MatrixIterator;

      VCT x( vec_ );  // Evaluation of the left-hand side sparse vector operand
      MCT A( mat_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( x.size()    == vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == mat_.columns(), "Invalid number of columns" );

      ElementType res = ElementType();

      VectorIterator velem( x.begin() );
      const VectorIterator vend( x.end() );
      if( velem == vend ) {
         reset( res );
         return res;
      }

      MatrixIterator melem( A.begin(index) );
      const MatrixIterator mend( A.end(index) );
      if( melem == mend ) {
         reset( res );
         return res;
      }

      while( true ) {
         if( velem->index() < melem->index() ) {
            ++velem;
            if( velem == vend ) break;
         }
         else if( melem->index() < velem->index() ) {
            ++melem;
            if( melem == mend ) break;
         }
         else {
            res = velem->value() * melem->value();
            ++velem;
            ++melem;
            break;
         }
      }

      if( melem != mend && velem != vend )
      {
         while( true ) {
            if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else {
               res += velem->value() * melem->value();
               ++velem;
               if( velem == vend ) break;
               ++melem;
               if( melem == mend ) break;
            }
         }
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

   //**NonZeros function***************************************************************************
   /*!\brief Returns an estimation for the number of non-zero elements in the sparse vector.
   //
   // \return The estimate for the number of non-zero elements in the sparse vector.
   */
   inline size_t nonZeros() const {
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
   /*!\brief Returns the right-hand side transpose sparse matrix operand.
   //
   // \return The right-hand side transpose sparse matrix operand.
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
      return ( vec_.isAliased( alias ) || mat_.isAliased( alias ) );
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
      return ( vec_.isAliased( alias ) || mat_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return ( size() > SMP_TSVECSMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vec_;  //!< Left-hand side sparse vector of the multiplication expression.
   RightOperand mat_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse vector-transpose sparse matrix multiplication to
   //        a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse vector-
   // transpose sparse matrix multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,true>& lhs, const TSVecTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Resetting the left-hand side target dense vector
      reset( ~lhs );

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the right-hand side sparse matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse vector-sparse matrix multiplication
      TSVecTSMatMultExpr::selectAssignKernel( ~lhs, x, A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Serial assignment to dense vectors**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serial assignment of a transpose sparse vector-transpose sparse matrix multiplication
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the serial assignment kernel for the transpose sparse vector-
   // transpose sparse matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<VT2>::Type::ConstIterator  VectorIterator;
      typedef typename RemoveReference<MT1>::Type::ConstIterator  MatrixIterator;

      const VectorIterator vend( x.end() );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         const MatrixIterator mend ( A.end(j)   );
         MatrixIterator       melem( A.begin(j) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         while( true ) {
            if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else {
               y[j] = velem->value() * melem->value();
               ++velem;
               ++melem;
               break;
            }
         }

         if( velem != vend && melem != mend )
         {
            while( true ) {
               if( velem->index() < melem->index() ) {
                  ++velem;
                  if( velem == vend ) break;
               }
               else if( melem->index() < velem->index() ) {
                  ++melem;
                  if( melem == mend ) break;
               }
               else {
                  y[j] += velem->value() * melem->value();
                  ++velem;
                  if( velem == vend ) break;
                  ++melem;
                  if( melem == mend ) break;
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose sparse vector-transpose sparse matrix multiplication
   //        (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the SMP assignment kernel for the transpose sparse vector-transpose
   // sparse matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      smpAssign( y, x * A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse vector-transpose sparse matrix multiplication to
   //        a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse vector-
   // transpose sparse matrix multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,true>& lhs, const TSVecTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename RemoveReference<LT>::Type::ConstIterator  VectorIterator;
      typedef typename RemoveReference<RT>::Type::ConstIterator  MatrixIterator;

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the right-hand side sparse matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse vector-sparse matrix multiplication
      ElementType accu;
      const VectorIterator vend( x.end() );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         const MatrixIterator mend ( A.end(j)   );
         MatrixIterator       melem( A.begin(j) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         reset( accu );

         while( true ) {
            if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else {
               accu = velem->value() * melem->value();
               ++velem;
               ++melem;
               break;
            }
         }

         if( velem != vend && melem != mend )
         {
            while( true ) {
               if( velem->index() < melem->index() ) {
                  ++velem;
                  if( velem == vend ) break;
               }
               else if( melem->index() < velem->index() ) {
                  ++melem;
                  if( melem == mend ) break;
               }
               else {
                  accu += velem->value() * melem->value();
                  ++velem;
                  if( velem == vend ) break;
                  ++melem;
                  if( melem == mend ) break;
               }
            }
         }

         if( !isDefault( accu ) )
            (~lhs).insert( j, accu );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose sparse vector-transpose sparse matrix multiplication
   //        to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose sparse
   // vector-transpose sparse matrix multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,true>& lhs, const TSVecTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the right-hand side sparse matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse matrix-sparse vector multiplication
      TSVecTSMatMultExpr::selectAddAssignKernel( ~lhs, x, A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Serial addition assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serial addition assignment of a transpose sparse vector-transpose sparse matrix
   //        multiplication (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the serial addition assignment kernel for the transpose sparse
   // vector-transpose sparse matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAddAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<VT2>::Type::ConstIterator  VectorIterator;
      typedef typename RemoveReference<MT1>::Type::ConstIterator  MatrixIterator;

      const VectorIterator vend( x.end() );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         const MatrixIterator mend ( A.end(j)   );
         MatrixIterator       melem( A.begin(j) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         while( true ) {
            if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else {
               y[j] += velem->value() * melem->value();
               ++velem;
               if( velem == vend ) break;
               ++melem;
               if( melem == mend ) break;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose sparse vector-transpose sparse matrix
   //        multiplication (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the SMP addition assignment kernel for the transpose sparse vector-
   // transpose sparse matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
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
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose sparse vector-transpose sparse matrix
   //        multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // sparse vector-transpose sparse matrix multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,true>& lhs, const TSVecTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the left-hand side sparse vector operand
      LT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the right-hand side sparse matrix operand
      RT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse matrix-sparse vector multiplication
      TSVecTSMatMultExpr::selectSubAssignKernel( ~lhs, x, A );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Serial subtraction assignment to dense vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serial subtraction assignment of a transpose sparse vector-transpose sparse matrix
   //        multiplication (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the serial subtraction assignment kernel for the transpose sparse
   // vector-transpose sparse matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectSubAssignKernel( VT1& y, const VT2& x, const MT1& A )
   {
      typedef typename RemoveReference<VT2>::Type::ConstIterator  VectorIterator;
      typedef typename RemoveReference<MT1>::Type::ConstIterator  MatrixIterator;

      const VectorIterator vend( x.end() );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         const MatrixIterator mend ( A.end(j)   );
         MatrixIterator       melem( A.begin(j) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         while( true ) {
            if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else {
               y[j] -= velem->value() * melem->value();
               ++velem;
               if( velem == vend ) break;
               ++melem;
               if( melem == mend ) break;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose sparse vector-transpose sparse matrix
   //        multiplication (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param x The left-hand side sparse vector operand.
   // \param A The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the SMP subtraction assignment kernel for the transpose sparse
   // vector-transpose sparse matrix multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename VT2    // Type of the left-hand side vector operand
           , typename MT1 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
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
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a transpose sparse vector-transpose sparse matrix
   //        multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // sparse vector-transpose sparse matrix multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,true>& lhs, const TSVecTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
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
//        column-major sparse matrix (\f$ \vec{y}^T=\vec{x}^T*A \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side transpose sparse vector for the multiplication.
// \param mat The right-hand side column-major sparse matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose sparse vector and a column-major
// sparse matrix:

   \code
   using blaze::rowVector;
   using blaze::columnMajor;

   blaze::CompressedVector<double,rowVector> x, y;
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// The operator returns an expression representing a transpose sparse vector of the higher-order
// element type of the two involved element types \a T1::ElementType and \a T2::ElementType.
// Both the sparse vector type \a T1 and the sparse matrix type \a T2 as well as the two element
// types \a T1::ElementType and \a T2::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of rows of
// the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side sparse vector
        , typename T2 >  // Type of the right-hand side sparse matrix
inline const typename DisableIf< IsMatMatMultExpr<T2>, TSVecTSMatMultExpr<T1,T2> >::Type
   operator*( const SparseVector<T1,true>& vec, const SparseMatrix<T2,true>& mat )
{
   BLAZE_FUNCTION_TRACE;

   if( (~vec).size() != (~mat).rows() )
      throw std::invalid_argument( "Vector and matrix sizes do not match" );

   return TSVecTSMatMultExpr<T1,T2>( ~vec, ~mat );
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
struct SubvectorExprTrait< TSVecTSMatMultExpr<VT,MT>, AF >
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
