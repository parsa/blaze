//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatSMatMultExpr.h
//  \brief Header file for the transpose dense matrix/sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATSMATMULTEXPR_H_


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
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseMatrix.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SMatDVecMultExprTrait.h>
#include <blaze/math/traits/SMatSVecMultExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/TDMatDVecMultExprTrait.h>
#include <blaze/math/traits/TDMatSVecMultExprTrait.h>
#include <blaze/math/traits/TDVecSMatMultExprTrait.h>
#include <blaze/math/traits/TDVecTDMatMultExprTrait.h>
#include <blaze/math/traits/TSVecTDMatMultExprTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-sparse matrix multiplications.
// \ingroup dense_matrix_expression
//
// The TDMatSMatMultExpr class represents the compile time expression for multiplications between
// a column-major dense matrix and a row-major sparse matrix.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
class TDMatSMatMultExpr : public DenseMatrix< TDMatSMatMultExpr<MT1,MT2>, true >
                        , private MatMatMultExpr
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT1::ResultType     RT1;  //!< Result type of the left-hand side dense matrix expression.
   typedef typename MT2::ResultType     RT2;  //!< Result type of the right-hand side sparse matrix expression.
   typedef typename RT1::ElementType    ET1;  //!< Element type of the left-hand side dense matrix expression.
   typedef typename RT2::ElementType    ET2;  //!< Element type of the right-hand side sparse matrix expression.
   typedef typename MT1::CompositeType  CT1;  //!< Composite type of the left-hand side dense matrix expression.
   typedef typename MT2::CompositeType  CT2;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   enum { evaluateLeft = IsComputation<MT1>::value || RequiresEvaluation<MT1>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse matrix expression.
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
   typedef TDMatSMatMultExpr<MT1,MT2>          This;           //!< Type of this TDMatSMatMultExpr instance.
   typedef typename MultTrait<RT1,RT2>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.
   typedef const ElementType                   ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                    CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT1>::value, const MT1, const MT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT2>::value, const MT2, const MT2& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side dense matrix operand.
   typedef typename SelectType< evaluateLeft, const RT1, CT1 >::Type  LT;

   //! Type for the assignment of the right-hand side sparse matrix operand.
   typedef typename SelectType< evaluateRight, const RT2, CT2 >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = !evaluateLeft && !evaluateRight };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDMatSMatMultExpr class.
   //
   // \param lhs The left-hand side dense matrix operand of the multiplication expression.
   // \param rhs The right-hand side sparse matrix operand of the multiplication expression.
   */
   explicit inline TDMatSMatMultExpr( const MT1& lhs, const MT2& rhs )
      : lhs_( lhs )  // Left-hand side dense matrix of the multiplication expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the multiplication expression
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

      if( lhs_.columns() != 0 ) {
         tmp = lhs_(i,0) * rhs_(0,j);
         for( size_t k=1; k<lhs_.columns(); ++k ) {
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
   /*!\brief Returns the left-hand side transpose dense matrix operand.
   //
   // \return The left-hand side transpose dense matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side sparse matrix operand.
   //
   // \return The right-hand side sparse matrix operand.
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
      return lhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return ( columns() > SMP_TDMATSMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-sparse matrix multiplication to a row-major
   //        dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,false>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectRowMajorAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to row-major dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a transpose dense matrix-sparse matrix
   // multiplication expression to a row-major dense matrix. This assign function is used in case
   // the element type of the target matrix is resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<MT3,MT4,MT5> >::Type
      selectRowMajorAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      for( size_t i=0; i<A.rows(); ++i ) {
         for( size_t j=0; j<C.columns(); ++j ) {
            reset( C(i,j) );
         }
         for( size_t j=0; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               if( isDefault( C(i,element->index()) ) )
                  C(i,element->index()) = A(i,j) * element->value();
               else
                  C(i,element->index()) += A(i,j) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to row-major dense matrices********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a row-major dense matrix. This assign function
   // is used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<MT3,MT4,MT5> >::Type
      selectRowMajorAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      const size_t last( A.rows() & size_t(-4) );

      for( size_t i=0; i<last; i+=4 ) {
         for( size_t j=0; j<C.columns(); ++j ) {
            reset( C(i  ,j) );
            reset( C(i+1,j) );
            reset( C(i+2,j) );
            reset( C(i+3,j) );
         }
         for( size_t j=0; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               C(i  ,element->index()) += A(i  ,j) * element->value();
               C(i+1,element->index()) += A(i+1,j) * element->value();
               C(i+2,element->index()) += A(i+2,j) * element->value();
               C(i+3,element->index()) += A(i+3,j) * element->value();
            }
         }
      }

      for( size_t i=last; i<A.rows(); ++i ) {
         for( size_t j=0; j<C.columns(); ++j ) {
            reset( C(i,j) );
         }
         for( size_t j=0; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               C(i,element->index()) += A(i,j) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to row-major dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP assignment of a transpose dense matrix-sparse matrix
   // multiplication expression to a row-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectRowMajorAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-sparse matrix multiplication to a column-major
   //        dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a sparse matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,true>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectColumnMajorAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to column-major dense matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a transpose dense matrix-sparse matrix
   // multiplication expression to a column-major dense matrix. This assign function is used in
   // case the element type of the target matrix is resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseDefaultKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      reset( C );

      for( size_t i=0; i<B.rows(); ++i ) {
         ConstIterator element( B.begin(i) );
         const ConstIterator end( B.end(i) );
         for( ; element!=end; ++element ) {
            for( size_t j=0; j<A.rows(); ++j ) {
               if( isDefault( C(j,element->index()) ) )
                  C(j,element->index()) = A(j,i) * element->value();
               else
                  C(j,element->index()) += A(j,i) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to column-major dense matrices*****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a column-major dense matrix. This assign function
   // is used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseOptimizedKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      reset( C );

      BLAZE_INTERNAL_ASSERT( ( A.rows() - ( A.rows() % 4UL ) ) == ( A.rows() & size_t(-4) ), "Invalid end calculation" );
      const size_t jend( A.rows() & size_t(-4) );
      size_t j( 0UL );

      for( size_t i=0UL; i<B.rows(); ++i ) {
         ConstIterator element( B.begin(i) );
         const ConstIterator end( B.end(i) );

         const size_t nonzeros( B.nonZeros(i) );
         const size_t kend( nonzeros & size_t(-4) );

         for( size_t k=0UL; k<kend; k+=4UL )
         {
            const size_t i1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t i2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t i3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t i4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            for( j=0UL; j<jend; j+=4UL ) {
               C(j    ,i1) += A(j    ,i) * v1;
               C(j+1UL,i1) += A(j+1UL,i) * v1;
               C(j+2UL,i1) += A(j+2UL,i) * v1;
               C(j+3UL,i1) += A(j+3UL,i) * v1;
               C(j    ,i2) += A(j    ,i) * v2;
               C(j+1UL,i2) += A(j+1UL,i) * v2;
               C(j+2UL,i2) += A(j+2UL,i) * v2;
               C(j+3UL,i2) += A(j+3UL,i) * v2;
               C(j    ,i3) += A(j    ,i) * v3;
               C(j+1UL,i3) += A(j+1UL,i) * v3;
               C(j+2UL,i3) += A(j+2UL,i) * v3;
               C(j+3UL,i3) += A(j+3UL,i) * v3;
               C(j    ,i4) += A(j    ,i) * v4;
               C(j+1UL,i4) += A(j+1UL,i) * v4;
               C(j+2UL,i4) += A(j+2UL,i) * v4;
               C(j+3UL,i4) += A(j+3UL,i) * v4;
            }
            for( ; j<A.rows(); ++j ) {
               C(j,i1) += A(j,i) * v1;
               C(j,i2) += A(j,i) * v2;
               C(j,i3) += A(j,i) * v3;
               C(j,i4) += A(j,i) * v4;
            }
         }

         for( ; element!=end; ++element ) {
            for( j=0UL; j<jend; j+=4UL ) {
               C(j    ,element->index()) += A(j    ,i) * element->value();
               C(j+1UL,element->index()) += A(j+1UL,i) * element->value();
               C(j+2UL,element->index()) += A(j+2UL,i) * element->value();
               C(j+3UL,element->index()) += A(j+3UL,i) * element->value();
            }
            for( ; j<A.rows(); ++j ) {
               C(j,element->index()) += A(j,i) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to column-major dense matrices***********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-sparse matrix multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP assignment of a transpose dense matrix-sparse matrix
   // multiplication expression to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-sparse matrix multiplication to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
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

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose dense matrix-sparse matrix multiplication to
   //        a row-major dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,false>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectRowMajorAddAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to row-major dense matrices***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectRowMajorAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      const size_t last( A.rows() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( A.rows() - ( A.rows() % 4UL ) ) == last, "Invalid end calculation" );

      for( size_t i=0UL; i<last; i+=4UL ) {
         for( size_t j=0UL; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               C(i    ,element->index()) += A(i    ,j) * element->value();
               C(i+1UL,element->index()) += A(i+1UL,j) * element->value();
               C(i+2UL,element->index()) += A(i+2UL,j) * element->value();
               C(i+3UL,element->index()) += A(i+3UL,j) * element->value();
            }
         }
      }

      for( size_t i=last; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               C(i,element->index()) += A(i,j) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to row-major dense matrices*****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP addition assignment of a transpose dense matrix-sparse
   // matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectRowMajorAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpAddAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to column-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose dense matrix-sparse matrix multiplication to
   //        a column-major dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,true>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectColumnMajorAddAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to column-major dense matrices********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      const size_t jend( A.rows() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( A.rows() - ( A.rows() % 4UL ) ) == jend, "Invalid end calculation" );
      size_t j( 0UL );

      for( size_t i=0UL; i<B.rows(); ++i ) {
         ConstIterator element( B.begin(i) );
         const ConstIterator end( B.end(i) );

         const size_t nonzeros( B.nonZeros(i) );
         const size_t kend( nonzeros & size_t(-4) );

         for( size_t k=0UL; k<kend; k+=4UL )
         {
            const size_t i1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t i2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t i3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t i4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            for( j=0UL; j<jend; j+=4UL ) {
               C(j    ,i1) += A(j    ,i) * v1;
               C(j+1UL,i1) += A(j+1UL,i) * v1;
               C(j+2UL,i1) += A(j+2UL,i) * v1;
               C(j+3UL,i1) += A(j+3UL,i) * v1;
               C(j    ,i2) += A(j    ,i) * v2;
               C(j+1UL,i2) += A(j+1UL,i) * v2;
               C(j+2UL,i2) += A(j+2UL,i) * v2;
               C(j+3UL,i2) += A(j+3UL,i) * v2;
               C(j    ,i3) += A(j    ,i) * v3;
               C(j+1UL,i3) += A(j+1UL,i) * v3;
               C(j+2UL,i3) += A(j+2UL,i) * v3;
               C(j+3UL,i3) += A(j+3UL,i) * v3;
               C(j    ,i4) += A(j    ,i) * v4;
               C(j+1UL,i4) += A(j+1UL,i) * v4;
               C(j+2UL,i4) += A(j+2UL,i) * v4;
               C(j+3UL,i4) += A(j+3UL,i) * v4;
            }
            for( ; j<A.rows(); ++j ) {
               C(j,i1) += A(j,i) * v1;
               C(j,i2) += A(j,i) * v2;
               C(j,i3) += A(j,i) * v3;
               C(j,i4) += A(j,i) * v4;
            }
         }

         for( ; element!=end; ++element ) {
            for( j=0UL; j<jend; j+=4UL ) {
               C(j    ,element->index()) += A(j    ,i) * element->value();
               C(j+1UL,element->index()) += A(j+1UL,i) * element->value();
               C(j+2UL,element->index()) += A(j+2UL,i) * element->value();
               C(j+3UL,element->index()) += A(j+3UL,i) * element->value();
            }
            for( ; j<A.rows(); ++j ) {
               C(j,element->index()) += A(j,i) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to column-major dense matrices**************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP addition assignment of a transpose dense matrix-sparse
   // matrix multiplication expression to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpAddAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to row-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-sparse matrix multiplication to a
   //        row-major dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,false>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectRowMajorSubAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to row-major dense matrices********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectRowMajorSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      const size_t last( A.rows() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( A.rows() - ( A.rows() % 4UL ) ) == last, "Invalid end calculation" );

      for( size_t i=0UL; i<last; i+=4UL ) {
         for( size_t j=0UL; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               C(i    ,element->index()) -= A(i    ,j) * element->value();
               C(i+1UL,element->index()) -= A(i+1UL,j) * element->value();
               C(i+2UL,element->index()) -= A(i+2UL,j) * element->value();
               C(i+3UL,element->index()) -= A(i+3UL,j) * element->value();
            }
         }
      }

      for( size_t i=last; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<B.rows(); ++j ) {
            ConstIterator element( B.begin(j) );
            const ConstIterator end( B.end(j) );
            for( ; element!=end; ++element ) {
               C(i,element->index()) -= A(i,j) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to row-major dense matrices**************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP subtraction assignment of a transpose dense matrix-sparse
   // matrix multiplication expression to a row-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectRowMajorSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      smpSubAssign( C, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to column-major dense matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-sparse matrix multiplication to a
   //        column-major dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,true>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectColumnMajorSubAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to column-major dense matrices*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename DisableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      typedef typename MT5::ConstIterator  ConstIterator;

      const size_t jend( A.rows() & size_t(-4) );
      BLAZE_INTERNAL_ASSERT( ( A.rows() - ( A.rows() % 4UL ) ) == jend, "Invalid end calculation" );
      size_t j( 0UL );

      for( size_t i=0UL; i<B.rows(); ++i ) {
         ConstIterator element( B.begin(i) );
         const ConstIterator end( B.end(i) );

         const size_t nonzeros( B.nonZeros(i) );
         const size_t kend( nonzeros & size_t(-4) );

         for( size_t k=0UL; k<kend; k+=4UL )
         {
            const size_t i1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t i2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t i3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t i4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            for( j=0UL; j<jend; j+=4UL ) {
               C(j    ,i1) -= A(j    ,i) * v1;
               C(j+1UL,i1) -= A(j+1UL,i) * v1;
               C(j+2UL,i1) -= A(j+2UL,i) * v1;
               C(j+3UL,i1) -= A(j+3UL,i) * v1;
               C(j    ,i2) -= A(j    ,i) * v2;
               C(j+1UL,i2) -= A(j+1UL,i) * v2;
               C(j+2UL,i2) -= A(j+2UL,i) * v2;
               C(j+3UL,i2) -= A(j+3UL,i) * v2;
               C(j    ,i3) -= A(j    ,i) * v3;
               C(j+1UL,i3) -= A(j+1UL,i) * v3;
               C(j+2UL,i3) -= A(j+2UL,i) * v3;
               C(j+3UL,i3) -= A(j+3UL,i) * v3;
               C(j    ,i4) -= A(j    ,i) * v4;
               C(j+1UL,i4) -= A(j+1UL,i) * v4;
               C(j+2UL,i4) -= A(j+2UL,i) * v4;
               C(j+3UL,i4) -= A(j+3UL,i) * v4;
            }
            for( ; j<A.rows(); ++j ) {
               C(j,i1) -= A(j,i) * v1;
               C(j,i2) -= A(j,i) * v2;
               C(j,i3) -= A(j,i) * v3;
               C(j,i4) -= A(j,i) * v4;
            }
         }

         for( ; element!=end; ++element ) {
            for( j=0UL; j<jend; j+=4UL ) {
               C(j    ,element->index()) -= A(j    ,i) * element->value();
               C(j+1UL,element->index()) -= A(j+1UL,i) * element->value();
               C(j+2UL,element->index()) -= A(j+2UL,i) * element->value();
               C(j+3UL,element->index()) -= A(j+3UL,i) * element->value();
            }
            for( ; j<A.rows(); ++j ) {
               C(j,element->index()) -= A(j,i) * element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to column-major dense matrices***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the SMP subtraction assignment of a transpose dense matrix-sparse
   // matrix multiplication expression to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline typename EnableIf< UseSMPAssignKernel<MT3,MT4,MT5> >::Type
      selectColumnMajorSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT2 );
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
/*!\brief Multiplication operator for the multiplication of a column-major dense matrix and a
//        row-major sparse matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a column-major dense matrix and a row-major
// sparse matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
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
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side sparse matrix
inline const TDMatSMatMultExpr<T1,T2>
   operator*( const DenseMatrix<T1,true>& lhs, const SparseMatrix<T2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   return TDMatSMatMultExpr<T1,T2>( ~lhs, ~rhs );
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
struct TDMatDVecMultExprTrait< TDMatSMatMultExpr<MT1,MT2>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsDenseVector<VT>::value   && IsColumnVector<VT>::value
                              , typename TDMatDVecMultExprTrait< MT1, typename SMatDVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename VT >
struct TDMatSVecMultExprTrait< TDMatSMatMultExpr<MT1,MT2>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsSparseVector<VT>::value  && IsColumnVector<VT>::value
                              , typename TDMatSVecMultExprTrait< MT1, typename SMatSVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT1, typename MT2 >
struct TDVecTDMatMultExprTrait< VT, TDMatSMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value   && IsRowVector<VT>::value          &&
                                IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value
                              , typename TDMatDVecMultExprTrait< MT1, typename SMatDVecMultExprTrait<MT2,VT>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT1, typename MT2 >
struct TSVecTDMatMultExprTrait< VT, TDMatSMatMultExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value  && IsRowVector<VT>::value          &&
                                IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value
                              , typename TDVecSMatMultExprTrait< typename TSVecTDMatMultExprTrait<VT,MT1>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool AF >
struct SubmatrixExprTrait< TDMatSMatMultExpr<MT1,MT2>, AF >
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
struct RowExprTrait< TDMatSMatMultExpr<MT1,MT2> >
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
struct ColumnExprTrait< TDMatSMatMultExpr<MT1,MT2> >
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
