//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecTDVecMultExpr.h
//  \brief Header file for the dense vector/dense vector outer product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECTDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECTDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECTDVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for outer products between two dense vectors.
// \ingroup dense_vector_expression
//
// The DVecTDVecMultExpr class represents the compile time expression for outer products
// between dense vectors.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
class DVecTDVecMultExpr : public DenseMatrix< DVecTDVecMultExpr<VT1,VT2>, false >
                        , private Expression
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT1::ResultType     RT1;  //!< Result type of the left-hand side dense vector expression.
   typedef typename VT2::ResultType     RT2;  //!< Result type of the right-hand side dense vector expression.
   typedef typename VT1::ReturnType     RN1;  //!< Return type of the left-hand side dense vector expression.
   typedef typename VT2::ReturnType     RN2;  //!< Return type of the right-hand side dense vector expression.
   typedef typename VT1::CompositeType  CT1;  //!< Composite type of the left-hand side dense vector expression.
   typedef typename VT2::CompositeType  CT2;  //!< Composite type of the right-hand side dense vector expression.
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
   typedef typename MultExprTrait<RN1,RN2>::Type  ExprReturnType;
   //**********************************************************************************************

   //**Evaluation strategy*************************************************************************
   //! Compilation switch for the evaluation strategy of the outer product expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch
       for the evaluation strategy of the outer product expression. In case either of the
       two dense vector operands is an expression, \a useAssign will be set to \a true and
       the outer product expression will be evaluated via the \a assign function family.
       Otherwise \a useAssign will be set to \a false and the expression will be evaluated
       via the function call operator. */
   enum { useAssign = ( IsComputation<VT1>::value || IsComputation<VT2>::value ) };

   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT >
   struct UseAssign {
      enum { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DVecTDVecMultExpr<VT1,VT2>          This;           //!< Type of this DVecTDVecMultExpr instance.
   typedef typename MultTrait<RT1,RT2>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DVecTDVecMultExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT1>::value, const VT1, const VT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT2>::value, const VT2, const VT2& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side dense vector operand.
   typedef typename SelectType< IsComputation<VT1>::value, const RT1, CT1 >::Type  LT;

   //! Type for the assignment of the right-hand side dense vector operand.
   typedef typename SelectType< IsComputation<VT2>::value, const RT2, CT2 >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecTDVecMultExpr class.
   //
   // \param lhs The left-hand side dense vector operand of the multiplication expression.
   // \param rhs The right-hand side dense vector operand of the multiplication expression.
   */
   explicit inline DVecTDVecMultExpr( const VT1& lhs, const VT2& rhs )
      : lhs_( lhs )  // Left-hand side dense vector of the multiplication expression
      , rhs_( rhs )  // Right-hand side dense vector of the multiplication expression
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < lhs_.size(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < rhs_.size(), "Invalid column access index" );

      return lhs_[i] * rhs_[j];
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const {
      return lhs_.size();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const {
      return rhs_.size();
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
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
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
      return ( lhs_.canAlias( alias ) || rhs_.canAlias( alias ) );
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
   LeftOperand  lhs_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side dense vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector outer product to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector outer product expression to a row-major dense matrix. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT> >::Type
      assign( DenseMatrix<MT,false>& lhs, const DVecTDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (~lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (~lhs).columns(), "Invalid vector size" );

      const size_t m( (~lhs).rows() );
      const size_t n( (~lhs).columns() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t jend( n & size_t(-2) );

      for( size_t i=0UL; i<m; ++i ) {
         for( size_t j=0UL; j<jend; j+=2UL ) {
            (~lhs)(i,j    ) = x[i] * y[j  ];
            (~lhs)(i,j+1UL) = x[i] * y[j+1];
         }
         if( jend < n ) {
            (~lhs)(i,jend) = x[i] * y[jend];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector outer product to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,true>& lhs, const DVecTDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (~lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (~lhs).columns(), "Invalid vector size" );

      const size_t m( (~lhs).rows() );
      const size_t n( (~lhs).columns() );

      BLAZE_INTERNAL_ASSERT( ( m - ( m % 2UL ) ) == ( m & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( m & size_t(-2) );

      for( size_t j=0UL; j<n; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            (~lhs)(i    ,j) = x[i  ] * y[j];
            (~lhs)(i+1UL,j) = x[i+1] * y[j];
         }
         if( iend < m ) {
            (~lhs)(iend,j) = x[iend] * y[j];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector outer product to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector outer product expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const DVecTDVecMultExpr& rhs )
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
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-dense vector outer product to a row-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // dense vector outer product expression to a row-major dense matrix. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT> >::Type
      addAssign( DenseMatrix<MT,false>& lhs, const DVecTDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (~lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (~lhs).columns(), "Invalid vector size" );

      const size_t m( (~lhs).rows() );
      const size_t n( (~lhs).columns() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t jend( n & size_t(-2) );

      for( size_t i=0UL; i<m; ++i ) {
         for( size_t j=0UL; j<jend; j+=2UL ) {
            (~lhs)(i,j    ) += x[i] * y[j    ];
            (~lhs)(i,j+1UL) += x[i] * y[j+1UL];
         }
         if( jend < n ) {
            (~lhs)(i,jend) += x[i] * y[jend];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to column-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-dense vector outer product to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // dense vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,true>& lhs, const DVecTDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (~lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (~lhs).columns(), "Invalid vector size" );

      const size_t m( (~lhs).rows() );
      const size_t n( (~lhs).columns() );

      BLAZE_INTERNAL_ASSERT( ( m - ( m % 2UL ) ) == ( m & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( m & size_t(-2) );

      for( size_t j=0UL; j<n; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            (~lhs)(i    ,j) += x[i    ] * y[j];
            (~lhs)(i+1UL,j) += x[i+1UL] * y[j];
         }
         if( iend < m ) {
            (~lhs)(iend,j) += x[iend] * y[j];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to row-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-dense vector outer product to a row-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // dense vector outer product expression to a row-major dense matrix. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT> >::Type
      subAssign( DenseMatrix<MT,false>& lhs, const DVecTDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (~lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (~lhs).columns(), "Invalid vector size" );

      const size_t m( (~lhs).rows() );
      const size_t n( (~lhs).columns() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t jend( n & size_t(-2) );

      for( size_t i=0UL; i<m; ++i ) {
         for( size_t j=0UL; j<jend; j+=2UL ) {
            (~lhs)(i,j    ) -= x[i] * y[j    ];
            (~lhs)(i,j+1UL) -= x[i] * y[j+1UL];
         }
         if( jend < n ) {
            (~lhs)(i,jend) -= x[i] * y[jend];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to column-major dense matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-dense vector outer product to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // dense vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,true>& lhs, const DVecTDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (~lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (~lhs).columns(), "Invalid vector size" );

      const size_t m( (~lhs).rows() );
      const size_t n( (~lhs).columns() );

      BLAZE_INTERNAL_ASSERT( ( m - ( m % 2UL ) ) == ( m & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( m & size_t(-2) );

      for( size_t j=0UL; j<n; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            (~lhs)(i    ,j) -= x[i    ] * y[j];
            (~lhs)(i+1UL,j) -= x[i+1UL] * y[j];
         }
         if( iend < m ) {
            (~lhs)(iend,j) -= x[iend] * y[j];
         }
      }
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE( VT2 );
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
/*!\brief Multiplication operator for the outer product of two dense vectors
//        (\f$ A=\vec{b}*\vec{c}^T \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector for the outer product.
// \param rhs The right-hand side transpose dense vector for the outer product.
// \return The resulting dense matrix.
//
// This operator represents the outer product between a dense vector and a transpose dense
// vector:

   \code
   using blaze::columnVector;
   using blaze::rowMajor;

   blaze::DynamicVector<double,columnVector> a, b;
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   A = a * trans(b);
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved element types \a T1::ElementType and \a T2::ElementType. Both
// dense vector types \a T1 and \a T2 as well as the two element types \a T1::ElementType
// and \a T2::ElementType have to be supported by the MultTrait class template.
*/
template< typename T1    // Type of the left-hand side dense vector
        , typename T2 >  // Type of the right-hand side dense vector
inline const DVecTDVecMultExpr<T1,T2>
   operator*( const DenseVector<T1,false>& lhs, const DenseVector<T2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return DVecTDVecMultExpr<T1,T2>( ~lhs, ~rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
