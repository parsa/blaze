//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatSMatAddExpr.h
//  \brief Header file for the dense matrix/sparse matrix addition expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSMATADDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSMATADDEXPR_H_


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
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/smp/DenseMatrix.h>
#include <blaze/math/smp/SparseMatrix.h>
#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATSMATADDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-sparse matrix additions.
// \ingroup dense_matrix_expression
//
// The DMatSMatAddExpr class represents the compile time expression for additions between
// a dense matrix and a sparse matrix with identical storage order.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , bool SO >     // Storage order
class DMatSMatAddExpr : public DenseMatrix< DMatSMatAddExpr<MT1,MT2,SO>, SO >
                      , private MatMatAddExpr
                      , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT1::ResultType  RT1;  //!< Result type of the left-hand side dense matrix expression.
   typedef typename MT2::ResultType  RT2;  //!< Result type of the right-hand side sparse matrix expression.
   typedef typename MT1::ReturnType  RN1;  //!< Return type of the left-hand side dense matrix expression.
   typedef typename MT2::ReturnType  RN2;  //!< Return type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either matrix operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   enum { returnExpr = !IsTemporary<RN1>::value && !IsTemporary<RN2>::value };

   //! Expression return type for the subscript operator.
   typedef typename AddExprTrait<RN1,RN2>::Type  ExprReturnType;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DMatSMatAddExpr<MT1,MT2,SO>         This;           //!< Type of this DMatSMatAddExpr instance.
   typedef typename AddTrait<RT1,RT2>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef const ResultType  CompositeType;

   //! Composite type of the left-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT1>::value, const MT1, const MT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT2>::value, const MT2, const MT2& >::Type  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatSMatAddExpr class.
   //
   // \param lhs The left-hand side dense matrix operand of the addition expression.
   // \param rhs The right-hand side sparse matrix operand of the addition expression.
   */
   explicit inline DMatSMatAddExpr( const MT1& lhs, const MT2& rhs )
      : lhs_( lhs )  // Left-hand side dense matrix of the addition expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the addition expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( lhs.columns() == rhs.columns(), "Invalid number of columns" );
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
      BLAZE_INTERNAL_ASSERT( j < lhs_.columns(), "Invalid column access index" );
      return lhs_(i,j) + rhs_(i,j);
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
      return lhs_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
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
      return ( IsExpression<MT1>::value && lhs_.canAlias( alias ) ) ||
             ( rhs_.canAlias( alias ) );
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

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the addition expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the addition expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-sparse matrix addition to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-sparse
   // matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO2>& lhs, const DMatSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsExpression<MT1>::value && (~lhs).isAliased( &rhs.lhs_ ) ) {
         smpAddAssign( ~lhs, rhs.rhs_ );
      }
      else {
         smpAssign   ( ~lhs, rhs.lhs_ );
         smpAddAssign( ~lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-sparse matrix addition to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-sparse
   // matrix addition expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO2>& lhs, const DMatSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      typedef typename SelectType< SO == SO2, ResultType, OppositeType >::Type  TmpType;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OppositeType, !SO );
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
   /*!\brief Addition assignment of a dense matrix-sparse matrix addition to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // sparse matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO2>& lhs, const DMatSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAddAssign( ~lhs, rhs.lhs_ );
      smpAddAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-sparse matrix addition to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix-
   // sparse matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO2>& lhs, const DMatSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSubAssign( ~lhs, rhs.lhs_ );
      smpSubAssign( ~lhs, rhs.rhs_ );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT1, MT2 );
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
/*!\brief Addition operator for the addition of a dense matrix and a sparse matrix
//        (\f$ A=B+C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the matrix addition.
// \param rhs The right-hand side sparse matrix to be added to the left-hand side matrix.
// \return The sum of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the addition of a dense matrix and a sparse matrix:

   \code
   blaze::DynamicMatrix<double> A, C;
   blaze::CompressedMatrix<double> B;
   // ... Resizing and initialization
   C = A + B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a T1::ElementType and \a T2::ElementType.
// Both matrix types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1  // Type of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side sparse matrix
        , bool SO >    // Storage order
inline const DMatSMatAddExpr<T1,T2,SO>
   operator+( const DenseMatrix<T1,SO>& lhs, const SparseMatrix<T2,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   return DMatSMatAddExpr<T1,T2,SO>( ~lhs, ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a sparse matrix and a dense matrix
//        (\f$ A=B+C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the matrix addition.
// \param rhs The right-hand side dense matrix to be added to the left-hand side matrix.
// \return The sum of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the addition of a sparse matrix and a dense matrix:

   \code
   blaze::CompressedMatrix<double> A;
   blaze::DynamicMatrix<double> B, C;
   // ... Resizing and initialization
   C = A + B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a T1::ElementType and \a T2::ElementType.
// Both matrix types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1  // Type of the left-hand side sparse matrix
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order
inline const DMatSMatAddExpr<T2,T1,SO>
   operator+( const SparseMatrix<T1,SO>& lhs, const DenseMatrix<T2,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   return DMatSMatAddExpr<T2,T1,SO>( ~rhs, ~lhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition operator for the addition of a dense matrix-sparse matrix addition
//        expression and a dense matrix (\f$ A=(B+C)+D \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-sparse matrix addition.
// \param rhs The right-hand side dense matrix.
// \return The sum of the two matrices.
//
// This operator implements a performance optimized treatment of the addition of a dense
// matrix-sparse matrix addition expression to a dense matrix.
*/
template< typename T1  // Type of the dense matrix of the left-hand side expression
        , typename T2  // Type of the sparse matrix of the left-hand side expression
        , bool SO1     // Storage order of the left-hand side expression
        , typename T3  // Type of the right-hand side dense matrix
        , bool SO2 >   // Storage order of the right-hand side dense matrix
inline const typename AddExprTrait< DMatSMatAddExpr<T1,T2,SO1>, T3 >::Type
   operator+( const DMatSMatAddExpr<T1,T2,SO1>& lhs, const DenseMatrix<T3,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() + (~rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction operator for the subtraction of a dense matrix-sparse matrix addition
//        expression and a dense matrix (\f$ A=(B+C)-D \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-sparse matrix addition.
// \param rhs The right-hand side dense matrix.
// \return The difference of the two matrices.
//
// This operator implements a performance optimized treatment of the subtraction of a dense
// matrix-sparse matrix addition expression and a dense matrix.
*/
template< typename T1  // Type of the dense matrix of the left-hand side expression
        , typename T2  // Type of the sparse matrix of the left-hand side expression
        , bool SO1     // Storage order of the left-hand side expression
        , typename T3  // Type of the right-hand side dense matrix
        , bool SO2 >   // Storage order of the right-hand side dense matrix
inline const typename SubExprTrait< DMatSMatAddExpr<T1,T2,SO1>, T3 >::Type
   operator-( const DMatSMatAddExpr<T1,T2,SO1>& lhs, const DenseMatrix<T3,SO2>& rhs )
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
template< typename MT1, typename MT2, typename MT3 >
struct DMatDMatAddExprTrait< DMatSMatAddExpr<MT1,MT2,false>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsDenseMatrix<MT3>::value  && IsRowMajorMatrix<MT3>::value
                              , typename DMatSMatAddExprTrait< typename DMatDMatAddExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct DMatTDMatAddExprTrait< DMatSMatAddExpr<MT1,MT2,false>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix <MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsDenseMatrix <MT3>::value && IsColumnMajorMatrix<MT3>::value
                              , typename DMatSMatAddExprTrait< typename DMatTDMatAddExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct TDMatDMatAddExprTrait< DMatSMatAddExpr<MT1,MT2,true>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix <MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsDenseMatrix <MT3>::value && IsRowMajorMatrix<MT3>::value
                              , typename DMatTSMatAddExprTrait< typename TDMatDMatAddExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct TDMatTDMatAddExprTrait< DMatSMatAddExpr<MT1,MT2,true>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix <MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsDenseMatrix <MT3>::value && IsColumnMajorMatrix<MT3>::value
                              , typename TDMatTSMatAddExprTrait< typename TDMatTDMatAddExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct DMatDMatSubExprTrait< DMatSMatAddExpr<MT1,MT2,false>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsDenseMatrix<MT3>::value  && IsRowMajorMatrix<MT3>::value
                              , typename DMatSMatAddExprTrait< typename DMatDMatSubExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct DMatTDMatSubExprTrait< DMatSMatAddExpr<MT1,MT2,false>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsDenseMatrix<MT3>::value  && IsColumnMajorMatrix<MT3>::value
                              , typename DMatSMatAddExprTrait< typename DMatTDMatSubExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct TDMatDMatSubExprTrait< DMatSMatAddExpr<MT1,MT2,true>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsDenseMatrix<MT3>::value  && IsRowMajorMatrix<MT3>::value
                              , typename DMatTSMatAddExprTrait< typename TDMatDMatSubExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct TDMatTDMatSubExprTrait< DMatSMatAddExpr<MT1,MT2,true>, MT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsDenseMatrix<MT3>::value  && IsColumnMajorMatrix<MT3>::value
                              , typename TDMatTSMatAddExprTrait< typename TDMatTDMatSubExprTrait<MT1,MT3>::Type, MT2 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SO, bool AF >
struct SubmatrixExprTrait< DMatSMatAddExpr<MT1,MT2,SO>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename AddExprTrait< typename SubmatrixExprTrait<const MT1,AF>::Type
                                , typename SubmatrixExprTrait<const MT2,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SO >
struct RowExprTrait< DMatSMatAddExpr<MT1,MT2,SO> >
{
 public:
   //**********************************************************************************************
   typedef typename AddExprTrait< typename RowExprTrait<const MT1>::Type
                                , typename RowExprTrait<const MT2>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SO >
struct ColumnExprTrait< DMatSMatAddExpr<MT1,MT2,SO> >
{
 public:
   //**********************************************************************************************
   typedef typename AddExprTrait< typename ColumnExprTrait<const MT1>::Type
                                , typename ColumnExprTrait<const MT2>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
