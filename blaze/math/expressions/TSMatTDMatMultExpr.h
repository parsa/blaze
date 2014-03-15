//=================================================================================================
/*!
//  \file blaze/math/expressions/TSMatTDMatMultExpr.h
//  \brief Header file for the transpose sparse matrix/transpose dense matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TSMATTDMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSMATTDMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseMatrix.h>
#include <blaze/math/smp/SparseMatrix.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/TDMatDVecMultExprTrait.h>
#include <blaze/math/traits/TDMatSVecMultExprTrait.h>
#include <blaze/math/traits/TSMatDVecMultExprTrait.h>
#include <blaze/math/traits/TDVecTDMatMultExprTrait.h>
#include <blaze/math/traits/TDVecTSMatMultExprTrait.h>
#include <blaze/math/traits/TSVecTDMatMultExprTrait.h>
#include <blaze/math/traits/TSVecTSMatMultExprTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATDMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose sparse matrix-transpose dense matrix multiplications.
// \ingroup dense_matrix_expression
//
// The TSMatTDMatMultExpr class represents the compile time expression for multiplications between
// a column-major sparse matrix and a column-major dense matrix.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
class TSMatTDMatMultExpr : public DenseMatrix< TSMatTDMatMultExpr<MT1,MT2>, true >
                         , private MatMatMultExpr
                         , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT1::ResultType     RT1;  //!< Result type of the left-hand side sparse matrix expression.
   typedef typename MT2::ResultType     RT2;  //!< Result type of the right-hand side dense matrix expression.
   typedef typename RT1::ElementType    ET1;  //!< Element type of the left-hand side dense matrix expression.
   typedef typename RT2::ElementType    ET2;  //!< Element type of the right-hand side sparse matrix expression.
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
   /*! In case the either of the two matrix operands requires an intermediate evaluation, the
       nested \value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseSMPAssignKernel {
      enum { value = evaluateLeft || evaluateRight };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case no SMP assignment is required and the element type of the target matrix has a
       fixed size (i.e. is not resizable), the nested \value will be set to 1, otherwise it
       will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseOptimizedKernel {
      enum { value = !UseSMPAssignKernel<T1,T2,T3>::value &&
                     !IsResizable<typename T1::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case no SMP assignment is required and the element type of the target matrix is
       resizable, the nested \value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseDefaultKernel {
      enum { value = !UseSMPAssignKernel<T1,T2,T3>::value &&
                     !UseOptimizedKernel<T1,T2,T3>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TSMatTDMatMultExpr<MT1,MT2>         This;           //!< Type of this TSMatTDMatMultExpr instance.
   typedef typename MultTrait<RT1,RT2>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.
   typedef const ElementType                   ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                    CompositeType;  //!< Data type for composite expression templates.

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
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateLeft && !evaluateRight };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TSMatTDMatMultExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the multiplication expression.
   // \param rhs The right-hand side dense matrix operand of the multiplication expression.
   */
   explicit inline TSMatTDMatMultExpr( const MT1& lhs, const MT2& rhs )
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

      ElementType tmp;

      if( lhs_.columns() != 0UL ) {
         tmp = lhs_(i,0UL) * rhs_(0UL,j);
         for( size_t k=1UL; k<lhs_.columns(); ++k ) {
            tmp += lhs_(i,k) * rhs_(k,j);
         }
      }
      else {
         reset( tmp );
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
   /*!\brief Returns the left-hand side transpose sparse matrix operand.
   //
   // \return The left-hand side transpose sparse matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side transpose dense matrix operand.
   //
   // \return The right-hand side transpose dense matrix operand.
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
      return ( columns() > SMP_TSMATTDMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-transpose dense matrix multiplication to a
   //        dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO>& lhs, const TSMatTDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the right-hand side sparse matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the left-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TSMatTDMatMultExpr::selectAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose sparse matrix-transpose dense matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a transpose sparse matrix-transpose
   // dense matrix multiplication expression to a dense matrix. This assign function is used
   // in case the element type of the target matrix is resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<MT3,MT4,MT5> >::Type
      selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      reset( C );

      const size_t block( 64UL );

      for( size_t jj=0UL; jj<B.columns(); jj+=block ) {
         const size_t jend( ( jj+block > B.columns() )?( B.columns() ):( jj+block ) );
         for( size_t i=0UL; i<A.columns(); ++i ) {
            const ConstIterator end( A.end(i) );
            ConstIterator element( A.begin(i) );
            for( ; element!=end; ++element ) {
               for( size_t j=jj; j<jend; ++j ) {
                  if( isDefault( C(element->index(),j) ) )
                     C(element->index(),j) = element->value() * B(i,j);
                  else
                     C(element->index(),j) += element->value() * B(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a transpose sparse matrix-transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // transpose dense matrix multiplication expression to a dense matrix. This assign function is
   // used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<MT3,MT4,MT5> >::Type
      selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      const size_t block( IsRowMajorMatrix<MT3>::value ? 128UL : 64UL );

      for( size_t jj=0UL; jj<B.columns(); jj+=block ) {
         const size_t jend( ( jj+block > B.columns() )?( B.columns() ):( jj+block ) );
         for( size_t i=0UL; i<A.rows(); ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               reset( C(i,j) );
            }
         }
         for( size_t i=0UL; i<A.columns(); ++i )
         {
            const ConstIterator end( A.end(i) );
            ConstIterator element( A.begin(i) );

            const size_t nonzeros( A.nonZeros(i) );

            const size_t kend( nonzeros & size_t(-4) );
            BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kend, "Invalid end calculation" );

            for( size_t k=0UL; k<kend; k+=4UL ) {
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

               for( size_t j=jj; j<jend; ++j ) {
                  C(i1,j) += v1 * B(i,j);
                  C(i2,j) += v2 * B(i,j);
                  C(i3,j) += v3 * B(i,j);
                  C(i4,j) += v4 * B(i,j);
               }
            }

            for( ; element!=end; ++element ) {
               for( size_t j=jj; j<jend; ++j ) {
                  C(element->index(),j) += element->value() * B(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose sparse matrix-transpose dense matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP assignment of a transpose sparse matrix-transpose dense
   // matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-transpose dense matrix multiplication to
   //        a sparse matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // transpose dense matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const TSMatTDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      typedef typename SelectType< SO, ResultType, OppositeType >::Type  TmpType;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename TmpType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( rhs );
      smpAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose sparse matrix-transpose dense matrix multiplication
   //        to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose sparse
   // matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO>& lhs, const TSMatTDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the right-hand side sparse matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the left-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TSMatTDMatMultExpr::selectAddAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a transpose sparse matrix-transpose dense matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose sparse
   // matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      const size_t block( IsRowMajorMatrix<MT3>::value ? 128UL : 64UL );

      for( size_t jj=0UL; jj<B.columns(); jj+=block ) {
         const size_t jend( ( jj+block > B.columns() )?( B.columns() ):( jj+block ) );
         for( size_t i=0UL; i<A.columns(); ++i )
         {
            const ConstIterator end( A.end(i) );
            ConstIterator element( A.begin(i) );

            const size_t nonzeros( A.nonZeros(i) );

            const size_t kend( nonzeros & size_t(-4) );
            BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kend, "Invalid end calculation" );

            for( size_t k=0UL; k<kend; k+=4UL ) {
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

               for( size_t j=jj; j<jend; ++j ) {
                  C(i1,j) += v1 * B(i,j);
                  C(i2,j) += v2 * B(i,j);
                  C(i3,j) += v3 * B(i,j);
                  C(i4,j) += v4 * B(i,j);
               }
            }

            for( ; element!=end; ++element ) {
               for( size_t j=jj; j<jend; ++j ) {
                  C(element->index(),j) += element->value() * B(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose sparse matrix-transpose dense matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP addition assignment of a transpose sparse matrix-transpose
   // dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpAddAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose sparse matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // sparse matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO>& lhs, const TSMatTDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the right-hand side sparse matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the left-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TSMatTDMatMultExpr::selectSubAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a transpose sparse matrix-transpose dense matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // sparse matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT4::ConstIterator  ConstIterator;

      const size_t block( IsRowMajorMatrix<MT3>::value ? 128UL : 64UL );

      for( size_t jj=0UL; jj<B.columns(); jj+=block ) {
         const size_t jend( ( jj+block > B.columns() )?( B.columns() ):( jj+block ) );
         for( size_t i=0UL; i<A.columns(); ++i )
         {
            const ConstIterator end( A.end(i) );
            ConstIterator element( A.begin(i) );

            const size_t nonzeros( A.nonZeros(i) );

            const size_t kend( nonzeros & size_t(-4) );
            BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kend, "Invalid end calculation" );

            for( size_t k=0UL; k<kend; k+=4UL ) {
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

               for( size_t j=jj; j<jend; ++j ) {
                  C(i1,j) -= v1 * B(i,j);
                  C(i2,j) -= v2 * B(i,j);
                  C(i3,j) -= v3 * B(i,j);
                  C(i4,j) -= v4 * B(i,j);
               }
            }

            for( ; element!=end; ++element ) {
               for( size_t j=jj; j<jend; ++j ) {
                  C(element->index(),j) -= element->value() * B(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose sparse matrix-transpose dense matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP subtraction assignment of a transpose sparse matrix-
   // transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpSubAssign( C, A * B );
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

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT2 );
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
/*!\brief Multiplication operator for the multiplication of a column-major sparse matrix and
//        a column-major dense matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the multiplication.
// \param rhs The right-hand side dense matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a column-major sparse matrix and a column-major
// dense matrix:

   \code
   using blaze::columnMajor;

   blaze::CompressedMatrix<double,columnMajor> A;
   blaze::DynamicMatrix<double,columnMajor> B, C;
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
inline const TSMatTDMatMultExpr<T1,T2>
   operator*( const SparseMatrix<T1,true>& lhs, const DenseMatrix<T2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   return TSMatTDMatMultExpr<T1,T2>( ~lhs, ~rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename VT >
struct TDMatDVecMultExprTrait< TSMatTDMatMultExpr<MT1,MT2>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value &&
                                IsDenseVector<VT>::value   && IsColumnVector<VT>::value
                              , typename TSMatDVecMultExprTrait< MT1, typename TDMatDVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename VT >
struct TDMatSVecMultExprTrait< TSMatTDMatMultExpr<MT1,MT2>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value &&
                                IsSparseVector<VT>::value  && IsColumnVector<VT>::value
                              , typename TSMatDVecMultExprTrait< MT1, typename TDMatSVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT1, typename MT2 >
struct TDVecTDMatMultExprTrait< VT, TSMatTDMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value   && IsRowVector<VT>::value          &&
                                IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value
                              , typename TDVecTDMatMultExprTrait< typename TDVecTSMatMultExprTrait<VT,MT1>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT1, typename MT2 >
struct TSVecTDMatMultExprTrait< VT, TSMatTDMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value  && IsRowVector<VT>::value          &&
                                IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value
                              , typename TSVecTDMatMultExprTrait< typename TSVecTSMatMultExprTrait<VT,MT1>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool AF >
struct SubmatrixExprTrait< TSMatTDMatMultExpr<MT1,MT2>, AF >
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
struct RowExprTrait< TSMatTDMatMultExpr<MT1,MT2> >
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
struct ColumnExprTrait< TSMatTDMatMultExpr<MT1,MT2> >
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
