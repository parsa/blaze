//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatSVecMultExpr.h
//  \brief Header file for the sparse matrix/sparse vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATSVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATSVECMULTEXPR_H_


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
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATDVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-sparse vector multiplications.
// \ingroup sparse_vector_expression
//
// The SMatSVecMultExpr class represents the compile time expression for multiplications
// between row-major sparse matrices and sparse vectors.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , typename VT >  // Type of the right-hand side sparse vector
class SMatSVecMultExpr : public SparseVector< SMatSVecMultExpr<MT,VT>, false >
                       , private MatVecMultExpr
                       , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     MRT;  //!< Result type of the left-hand side sparse matrix expression.
   typedef typename VT::ResultType     VRT;  //!< Result type of the right-hand side sparse vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the left-hand side sparse matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the right-hand side sparse vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side sparse matrix expression.
   enum { evaluateMatrix = RequiresEvaluation<MT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse vector expression.
   enum { evaluateVector = IsComputation<VT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the either the matrix or the vector operand require an intermediate evaluation,
       the nested \value will be set to 0, otherwise it will be 1. */
   template< typename T1, typename T2, typename T3 >
   struct UseSMPAssignKernel {
      enum { value = evaluateMatrix || evaluateVector };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SMatSVecMultExpr<MT,VT>             This;           //!< Type of this SMatSVecMultExpr instance.
   typedef typename MultTrait<MRT,VRT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.
   typedef const ElementType                   ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef const ResultType  CompositeType;

   //! Composite type of the left-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side sparse matrix operand.
   typedef MCT  LT;

   //! Type for the assignment of the right-hand side sparse vector operand.
   typedef typename SelectType< evaluateVector, const VRT, VCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateMatrix && !evaluateVector };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatSVecMultExpr class.
   //
   // \param mat The left-hand side sparse matrix operand of the multiplication expression.
   // \param vec The right-hand side sparse vector operand of the multiplication expression.
   */
   explicit inline SMatSVecMultExpr( const MT& mat, const VT& vec )
      : mat_( mat )  // Left-hand side sparse matrix of the multiplication expression
      , vec_( vec )  // Right-hand side sparse vector of the multiplication expression
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

      typedef typename RemoveReference<MCT>::Type::ConstIterator  MatrixIterator;
      typedef typename RemoveReference<VCT>::Type::ConstIterator  VectorIterator;

      ElementType res = ElementType();

      // Early exit
      if( vec_.size() == 0UL )
         return res;

      // Fast computation in case both the left-hand side matrix operand and the right-hand
      // side vector operand directly provide iterators
      if( !RequiresEvaluation<MT>::value && !RequiresEvaluation<VT>::value )
      {
         MCT A( mat_ );  // Evaluation of the left-hand side sparse matrix operand
         VCT x( vec_ );  // Evaluation of the right-hand side sparse vector operand

         MatrixIterator melem( A.begin(index) );
         const MatrixIterator mend( A.end(index) );
         if( melem == mend ) {
            return res;
         }

         VectorIterator velem( x.begin() );
         const VectorIterator vend( x.end() );
         if( velem == vend ) {
            return res;
         }

         while( true ) {
            if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else {
               res = melem->value() * velem->value();
               ++melem;
               ++velem;
               break;
            }
         }

         if( melem != mend && velem != vend )
         {
            while( true ) {
               if( melem->index() < velem->index() ) {
                  ++melem;
                  if( melem == mend ) break;
               }
               else if( velem->index() < melem->index() ) {
                  ++velem;
                  if( velem == vend ) break;
               }
               else {
                  res += melem->value() * velem->value();
                  ++melem;
                  if( melem == mend ) break;
                  ++velem;
                  if( velem == vend ) break;
               }
            }
         }
      }

      // Optimized computation in case the left-hand side matrix operand directly provides iterators
      else if( !RequiresEvaluation<MT>::value )
      {
         MCT A( mat_ );  // Evaluation of the left-hand side sparse matrix operand

         MatrixIterator melem( A.begin(index) );
         const MatrixIterator mend( A.end(index) );

         if( melem == mend )
            return res;

         res = melem->value() * vec_[melem->index()];
         ++melem;
         for( ; melem!=mend; ++melem ) {
            res += melem->value() * vec_[melem->index()];
         }
      }

      // Optimized computation in case the right-hand side vector operand directly provides iterators
      else if( !RequiresEvaluation<VT>::value )
      {
         VCT x( vec_ );  // Evaluation of the right-hand side sparse vector operand

         VectorIterator velem( x.begin() );
         const VectorIterator vend( x.end() );

         if( velem == vend )
            return res;

         res = mat_(index,velem->index()) * velem->value();
         ++velem;
         for( ; velem!=vend; ++velem ) {
            res += mat_(index,velem->index()) * velem->value();
         }
      }

      // Default computation in case both operands don't provide iterators
      else {
         res = mat_(index,0UL) * vec_[0UL];
         for( size_t j=1UL; j<vec_.size(); ++j ) {
            res += mat_(index,j) * vec_[j];
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
      return mat_.rows();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns an estimation for the number of non-zero elements in the sparse vector.
   //
   // \return The estimate for the number of non-zero elements in the sparse vector.
   */
   inline size_t nonZeros() const {
      return mat_.rows();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse matrix operand.
   //
   // \return The left-hand side sparse matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return mat_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side sparse vector operand.
   //
   // \return The right-hand side sparse vector operand.
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
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return ( size() > SMP_SMATSVECMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side sparse matrix of the multiplication expression.
   RightOperand vec_;  //!< Right-hand side sparse vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-sparse
   // vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const SMatSVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Resetting the left-hand side target dense vector
      reset( ~lhs );

      // Evaluation of the right-hand side sparse vector operand
      RT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side sparse matrix operand
      LT A( rhs.mat_ );

      // Checking the evaluated operators
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse matrix-sparse vector multiplication
      SMatSVecMultExpr::selectAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Serial assignment to dense vectors**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serial assignment of a sparse matrix-sparse vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side sparse matrix operand.
   // \param x The right-hand side sparse vector operand.
   // \return void
   //
   // This function implements the serial assignment kernel for the sparse matrix-sparse vector
   // multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  MatrixIterator;
      typedef typename RemoveReference<RT>::Type::ConstIterator  VectorIterator;

      const VectorIterator vend( x.end() );

      for( size_t i=0UL; i<y.size(); ++i )
      {
         const MatrixIterator mend ( A.end(i)   );
         MatrixIterator       melem( A.begin(i) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         while( true ) {
            if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else {
               y[i] = melem->value() * velem->value();
               ++melem;
               ++velem;
               break;
            }
         }

         if( melem != mend && velem != vend )
         {
            while( true ) {
               if( melem->index() < velem->index() ) {
                  ++melem;
                  if( melem == mend ) break;
               }
               else if( velem->index() < melem->index() ) {
                  ++velem;
                  if( velem == vend ) break;
               }
               else {
                  y[i] += melem->value() * velem->value();
                  ++melem;
                  if( melem == mend ) break;
                  ++velem;
                  if( velem == vend ) break;
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse matrix-sparse vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side sparse matrix operand.
   // \param x The right-hand side sparse vector operand.
   // \return void
   //
   // This function implements the SMP assignment kernel for the sparse matrix-sparse vector
   // multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      smpAssign( y, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-sparse vector multiplication to a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-sparse
   // vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const SMatSVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename RemoveReference<LT>::Type::ConstIterator  MatrixIterator;
      typedef typename RemoveReference<RT>::Type::ConstIterator  VectorIterator;

      RT x( rhs.vec_ );  // Evaluation of the right-hand side sparse vector operand
      if( x.nonZeros() == 0UL ) return;

      LT A( rhs.mat_ );  // Evaluation of the left-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      ElementType accu;
      const VectorIterator vend( x.end() );

      for( size_t i=0UL; i<(~lhs).size(); ++i )
      {
         const MatrixIterator mend ( A.end(i)   );
         MatrixIterator       melem( A.begin(i) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         reset( accu );

         while( true ) {
            if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else {
               accu = melem->value() * velem->value();
               ++melem;
               ++velem;
               break;
            }
         }

         if( melem != mend && velem != vend )
         {
            while( true ) {
               if( melem->index() < velem->index() ) {
                  ++melem;
                  if( melem == mend ) break;
               }
               else if( velem->index() < melem->index() ) {
                  ++velem;
                  if( velem == vend ) break;
               }
               else {
                  accu += melem->value() * velem->value();
                  ++melem;
                  if( melem == mend ) break;
                  ++velem;
                  if( velem == vend ) break;
               }
            }
         }

         if( !isDefault( accu ) )
            (~lhs).insert( i, accu );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // sparse vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const SMatSVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side sparse matrix operand
      LT A( rhs.mat_ );

      // Checking the evaluated operators
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse matrix-sparse vector multiplication
      SMatSVecMultExpr::selectAddAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Serial addition assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serial addition assignment of a sparse matrix-sparse vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side sparse matrix operand.
   // \param x The right-hand side sparse vector operand.
   // \return void
   //
   // This function implements the serial addition assignment kernel for the sparse matrix-sparse
   // vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  MatrixIterator;
      typedef typename RemoveReference<RT>::Type::ConstIterator  VectorIterator;

      const VectorIterator vend( x.end() );

      for( size_t i=0UL; i<y.size(); ++i )
      {
         const MatrixIterator mend ( A.end(i)   );
         MatrixIterator       melem( A.begin(i) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         while( true ) {
            if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else {
               y[i] += melem->value() * velem->value();
               ++melem;
               if( melem == mend ) break;
               ++velem;
               if( velem == vend ) break;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix-sparse vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side sparse matrix operand.
   // \param x The right-hand side sparse vector operand.
   // \return void
   //
   // This function implements the SMP addition assignment kernel for the sparse matrix-sparse
   // vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      smpAddAssign( y, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // sparse vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const SMatSVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename RemoveReference<LT>::Type::ConstIterator  MatrixIterator;
      typedef typename RemoveReference<RT>::Type::ConstIterator  VectorIterator;

      // Evaluation of the right-hand side sparse vector operand
      RT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side sparse matrix operand
      LT A( rhs.mat_ );

      // Checking the evaluated operators
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      // Performing the sparse matrix-sparse vector multiplication
      SMatSVecMultExpr::selectSubAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Serial subtraction assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serial subtraction assignment of a sparse matrix-sparse vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side sparse matrix operand.
   // \param x The right-hand side sparse vector operand.
   // \return void
   //
   // This function implements the serial subtraction assignment kernel for the sparse matrix-
   // sparse vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef typename RemoveReference<LT>::Type::ConstIterator  MatrixIterator;
      typedef typename RemoveReference<RT>::Type::ConstIterator  VectorIterator;

      const VectorIterator vend( x.end() );

      for( size_t i=0UL; i<y.size(); ++i )
      {
         const MatrixIterator mend ( A.end(i)   );
         MatrixIterator       melem( A.begin(i) );

         if( melem == mend ) continue;

         VectorIterator velem( x.begin() );

         while( true ) {
            if( melem->index() < velem->index() ) {
               ++melem;
               if( melem == mend ) break;
            }
            else if( velem->index() < melem->index() ) {
               ++velem;
               if( velem == vend ) break;
            }
            else {
               y[i] -= melem->value() * velem->value();
               ++melem;
               if( melem == mend ) break;
               ++velem;
               if( velem == vend ) break;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix-sparse vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup sparse_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side sparse matrix operand.
   // \param x The right-hand side sparse vector operand.
   // \return void
   //
   // This function implements the SMP subtraction assignment kernel for the sparse matrix-sparse
   // vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSMPAssignKernel<VT1,MT1,VT2> >::Type
      selectSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      smpSubAssign( y, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse matrix-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // matrix-sparse vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const SMatSVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT );
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
//        sparse vector (\f$ \vec{a}=B*\vec{c} \f$).
// \ingroup sparse_vector
//
// \param mat The left-hand side sparse matrix for the multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a row-major sparse matrix and a sparse
// vector:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::CompressedMatrix<double,rowMajor> A;
   blaze::CompressedVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns an expression representing a sparse vector of the higher-order element
// type of the two involved element types \a T1::ElementType and \a T2::ElementType. Both the
// sparse matrix type \a T1 and the sparse vector type \a T2 as well as the two element types
// \a T1::ElementType and \a T2::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side sparse vector
inline const typename DisableIf< IsMatMatMultExpr<T1>, SMatSVecMultExpr<T1,T2> >::Type
   operator*( const SparseMatrix<T1,false>& mat, const SparseVector<T2,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   if( (~mat).columns() != (~vec).size() )
      throw std::invalid_argument( "Matrix and vector sizes do not match" );

   return SMatSVecMultExpr<T1,T2>( ~mat, ~vec );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a sparse matrix-matrix
//        multiplication expression and a sparse vector (\f$ \vec{y}=(A*B)*\vec{x} \f$).
// \ingroup sparse_vector
//
// \param mat The left-hand side sparse matrix-matrix multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a sparse
// matrix-matrix multiplication expression and a sparse vector. It restructures the expression
// \f$ \vec{y}=(A*B)*\vec{x} \f$ to the expression \f$ \vec{y}=A*(B*\vec{x}) \f$.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side sparse vector
inline const typename EnableIf< IsMatMatMultExpr<T1>, MultExprTrait<T1,T2> >::Type::Type
   operator*( const SparseMatrix<T1,SO>& mat, const SparseVector<T2,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return (~mat).leftOperand() * ( (~mat).rightOperand() * vec );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT, bool AF >
struct SubvectorExprTrait< SMatSVecMultExpr<MT,VT>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename MultExprTrait< typename SubmatrixExprTrait<const MT,AF>::Type, VT >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
