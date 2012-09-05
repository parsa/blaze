//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatScalarMultExpr.h
//  \brief Header file for the dense matrix/scalar multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSCALARMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSCALARMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Expression.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATSCALARMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-scalar multiplications.
// \ingroup dense_matrix_expression
//
// The DMatScalarMultExpr class represents the compile time expression for multiplications between
// a dense matrix and a scalar value.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar value
        , bool SO >    // Storage order
class DMatScalarMultExpr : public DenseMatrix< DMatScalarMultExpr<MT,ST,SO>, SO >
                         , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     RT;  //!< Result type of the dense matrix expression.
   typedef typename MT::ReturnType     RN;  //!< Return type of the dense matrix expression.
   typedef typename MT::CompositeType  CT;  //!< Composite type of the dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case the dense matrix
       operand requires an intermediate evaluation, \a useAssign will be set to \a true and
       the multiplication expression will be evaluated via the \a assign function family.
       Otherwise \a useAssign will be set to \a false and the expression will be evaluated
       via the subscript operator. */
   enum { useAssign = !IsReference<CT>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct UseAssign {
      enum { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DMatScalarMultExpr<MT,ST,SO>               This;           //!< Type of this DMatScalarMultExpr instance.
   typedef typename MathTrait<RT,ST>::MultType        ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType          OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType         TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType           ElementType;    //!< Resulting element type.
   typedef const typename MultExprTrait<RN,ST>::Type  ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DMatScalarMultExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  LeftOperand;

   //! Composite type of the right-hand side scalar value.
   typedef typename MathTrait< typename BaseElementType<MT>::Type, ST >::MultType  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<MT>::value };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatScalarMultExpr class.
   //
   // \param matrix The left-hand side dense matrix of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   explicit inline DMatScalarMultExpr( const MT& matrix, ST scalar )
      : matrix_( matrix )  // Left-hand side dense matrix of the multiplication expression
      , scalar_( scalar )  // Right-hand side scalar of the multiplication expression
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
      BLAZE_INTERNAL_ASSERT( i < matrix_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < matrix_.columns(), "Invalid column access index" );
      return matrix_(i,j) * scalar_;
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const {
      return matrix_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const {
      return matrix_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return matrix_;
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
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return matrix_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  matrix_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a row-major dense matrix-scalar multiplication to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // multiplication expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( DenseMatrix<MT2,false>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );

      const size_t m( (~lhs).rows()    );
      const size_t n( (~lhs).columns() );
      for( size_t i=0UL; i<m; ++i ) {
         for( size_t j=0UL; j<n; ++j ) {
            (~lhs)(i,j) *= rhs.scalar_;
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a column-major dense matrix-scalar multiplication to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // multiplication expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( DenseMatrix<MT2,true>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );

      const size_t m( (~lhs).rows()    );
      const size_t n( (~lhs).columns() );
      for( size_t j=0UL; j<n; ++j ) {
         for( size_t i=0UL; i<m; ++i ) {
            (~lhs)(i,j) *= rhs.scalar_;
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a row-major dense matrix-scalar multiplication to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // multiplication expression to a sparse matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target sparse matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( SparseMatrix<MT2,false>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );

      for( size_t i=0UL; i<(~lhs).rows(); ++i )
      {
         typename MT2::Iterator element( (~lhs).begin(i) );
         const typename MT2::Iterator end( (~lhs).end(i) );

         for( ; element!=end; ++element )
            element->value() *= rhs.scalar_;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a column-major dense matrix-scalar multiplication to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // multiplication expression to a sparse matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target sparse matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( SparseMatrix<MT2,true>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );

      for( size_t j=0UL; j<(~lhs).columns(); ++j )
      {
         typename MT2::Iterator element( (~lhs).begin(j) );
         const typename MT2::Iterator end( (~lhs).end(j) );

         for( ; element!=end; ++element )
            element->value() *= rhs.scalar_;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-scalar multiplication to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // scalar multiplication expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case the
   // matrix operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      addAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-scalar multiplication to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix-
   // scalar multiplication expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      subAssign( ~lhs, tmp );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL UNARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Unary minus operator for the negation of a dense matrix (\f$ A = -B \f$).
// \ingroup dense_matrix
//
// \param dm The dense matrix to be negated.
// \return The negation of the matrix.
//
// This operator represents the negation of a dense matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = -A;
   \endcode

// The operator returns an expression representing the negation of the given dense matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatScalarMultExpr<MT,typename BaseElementType<MT>::Type,SO>
   operator-( const DenseMatrix<MT,SO>& dm )
{
   typedef typename BaseElementType<MT>::Type  ElementType;
   return DMatScalarMultExpr<MT,ElementType,SO>( ~dm, ElementType(-1) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense matrix and a scalar value
//        (\f$ A=B*s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result matrix.
//
// This operator represents the multiplication between a dense matrix and a scalar value:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = A * 1.25;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a T1::ElementType and \a T2. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename T1    // Type of the left-hand side dense matrix
        , bool SO        // Storage order of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsNumeric<T2>, typename MultExprTrait<T1,T2>::Type >::Type
   operator*( const DenseMatrix<T1,SO>& mat, T2 scalar )
{
   typedef typename MultExprTrait<T1,T2>::Type  Type;
   return Type( ~mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense matrix
//        (\f$ A=s*B \f$).
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param mat The right-hand side dense matrix for the multiplication.
// \return The scaled result matrix.
//
// This operator represents the multiplication between a a scalar value and dense matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = 1.25 * A;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a T1 and \a T2::ElementType. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order of the right-hand side dense matrix
inline const typename EnableIf< IsNumeric<T1>, typename MultExprTrait<T1,T2>::Type >::Type
   operator*( T1 scalar, const DenseMatrix<T2,SO>& mat )
{
   typedef typename MultExprTrait<T1,T2>::Type  Type;
   return Type( ~mat, scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar multiplication
//        expression and a scalar value (\f$ A=(B*s1)*s2 \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the multiplication of a
// dense matrix-scalar multiplication expression and a scalar value.
*/
template< typename MT     // Type of the dense matrix of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool SO         // Storage order of the dense matrix
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsNumeric<ST2>
                              , typename MultExprTrait< DMatScalarMultExpr<MT,ST1,SO>, ST2 >::Type >::Type
   operator*( const DMatScalarMultExpr<MT,ST1,SO>& mat, ST2 scalar )
{
   return mat.leftOperand() * ( mat.rightOperand() * scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense matrix-
//        scalar multiplication expression (\f$ A=s2*(B*s1) \f$).
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param mat The right-hand side dense matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a dense matrix-scalar multiplication expression.
*/
template< typename ST1  // Type of the left-hand side scalar
        , typename MT   // Type of the dense matrix of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO >     // Storage order of the dense matrix
inline const typename EnableIf< IsNumeric<ST1>
                              , typename MultExprTrait< ST1, DMatScalarMultExpr<MT,ST2,SO> >::Type >::Type
   operator*( ST1 scalar, const DMatScalarMultExpr<MT,ST2,SO>& mat )
{
   return mat.leftOperand() * ( scalar * mat.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a dense matrix-scalar multiplication
//        expression by a scalar value (\f$ A=(B*s1)/s2 \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the division of a
// dense matrix-scalar multiplication expression by a scalar value.
*/
template< typename MT     // Type of the dense matrix of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool SO         // Storage order of the dense matrix
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsFloatingPoint<typename MathTrait<ST1,ST2>::DivType>
                              , typename DivExprTrait< DMatScalarMultExpr<MT,ST1,SO>, ST2 >::Type >::Type
   operator/( const DMatScalarMultExpr<MT,ST1,SO>& mat, ST2 scalar )
{
   return mat.leftOperand() * ( mat.rightOperand() / scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar
//        multiplication expression and a dense vector (\f$ \vec{a}=(B*s1)*\vec{c} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param vec The right-hand side dense vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense matrix-scalar multiplication and a dense vector. It restructures the expression
// \f$ \vec{a}=(B*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Type of the dense matrix of the left-hand side expression
        , typename ST    // Type of the scalar of the left-hand side expression
        , bool SO        // Storage order of the left-hand side expression
        , typename VT >  // Type of the right-hand side dense vector
inline const typename MultExprTrait< DMatScalarMultExpr<MT,ST,SO>, VT >::Type
   operator*( const DMatScalarMultExpr<MT,ST,SO>& mat, const DenseVector<VT,false>& vec )
{
   return ( mat.leftOperand() * (~vec) ) * mat.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector and a dense
//        matrix-scalar multiplication expression (\f$ \vec{a}^T=\vec{c}^T*(B*s1) \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector.
// \param mat The right-hand side dense matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense vector and a dense matrix-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{c}^T*(B*s1) \f$ to the expression \f$ \vec{a}^T=(\vec{c}^T*B)*s1 \f$.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename MT  // Type of the dense matrix of the right-hand side expression
        , typename ST  // Type of the scalar of the right-hand side expression
        , bool SO >    // Storage order of the right-hand side expression
inline const typename MultExprTrait< VT, DMatScalarMultExpr<MT,ST,SO> >::Type
   operator*( const DenseVector<VT,true>& vec, const DMatScalarMultExpr<MT,ST,SO>& mat )
{
   return ( (~vec) * mat.leftOperand() ) * mat.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar
//        multiplication expression and a dense vector-scalar multiplication expression
//        (\f$ \vec{a}=(B*s1)*(\vec{c}*s2) \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param vec The right-hand side dense vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense matrix-scalar multiplication and a dense vector-scalar multiplication. It
// restructures the expression \f$ \vec{a}=(B*s1)*(\vec{c}*s2) \f$ to the expression
// \f$ \vec{a}=(B*\vec{c})*(s1*s2) \f$.
*/
template< typename MT     // Type of the dense matrix of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool SO         // Storage order of the left-hand side expression
        , typename VT     // Type of the dense vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar of the right-hand side expression
inline const DVecScalarMultExpr<typename MultExprTrait<MT,VT>::Type,typename MathTrait<ST1,ST2>::MultType,false>
   operator*( const DMatScalarMultExpr<MT,ST1,SO>& mat, const DVecScalarMultExpr<VT,ST2,false>& vec )
{
   return ( mat.leftOperand() * vec.leftOperand() ) * ( mat.rightOperand() * vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector-scalar
//        multiplication expression and a dense matrix-scalar multiplication expression
//        (\f$ \vec{a}^T=\vec{b}^T*(C*s1) \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector-scalar multiplication.
// \param mat The right-hand side dense matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense vector-scalar multiplication and a dense matrix-scalar multiplication. It
// restructures the expression \f$ \vec{a}=(\vec{b}^T*s1)*(C*s2) \f$ to the expression
// \f$ \vec{a}^T=(\vec{b}^T*C)*(s1*s2) \f$.
*/
template< typename VT   // Type of the dense vector of the left-hand side expression
        , typename ST1  // Type of the scalar of the left-hand side expression
        , typename MT   // Type of the dense matrix of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO >     // Storage order of the right-hand side expression
inline const typename MultExprTrait< DVecScalarMultExpr<VT,ST1,true>, DMatScalarMultExpr<MT,ST2,SO> >::Type
   operator*( const DVecScalarMultExpr<VT,ST1,true>& vec, const DMatScalarMultExpr<MT,ST2,SO>& mat )
{
   return ( vec.leftOperand() * mat.leftOperand() ) * ( vec.rightOperand() * mat.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar
//        multiplication expression and a sparse vector (\f$ \vec{a}=(B*s1)*\vec{c} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param vec The right-hand side sparse vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense matrix-scalar multiplication and a sparse vector. It restructures the expression
// \f$ \vec{a}=(B*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Type of the dense matrix of the left-hand side expression
        , typename ST    // Type of the scalar of the left-hand side expression
        , bool SO        // Storage order of the left-hand side expression
        , typename VT >  // Type of the right-hand side sparse vector
inline const typename MultExprTrait< DMatScalarMultExpr<MT,ST,SO>, VT >::Type
   operator*( const DMatScalarMultExpr<MT,ST,SO>& mat, const SparseVector<VT,false>& vec )
{
   return ( mat.leftOperand() * (~vec) ) * mat.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse vector and a dense
//        matrix-scalar multiplication expression (\f$ \vec{a}^T=\vec{c}^T*(B*s1) \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side sparse vector.
// \param mat The right-hand side dense matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a
// sparse vector and a dense matrix-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{c}^T*(B*s1) \f$ to the expression \f$ \vec{a}^T=(\vec{c}^T*B)*s1 \f$.
*/
template< typename VT  // Type of the left-hand side sparse vector
        , typename MT  // Type of the dense matrix of the right-hand side expression
        , typename ST  // Type of the scalar of the right-hand side expression
        , bool SO >    // Storage order of the right-hand side expression
inline const typename MultExprTrait< VT, DMatScalarMultExpr<MT,ST,SO> >::Type
   operator*( const SparseVector<VT,true>& vec, const DMatScalarMultExpr<MT,ST,SO>& mat )
{
   return ( (~vec) * mat.leftOperand() ) * mat.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar
//        multiplication expression and a sparse vector-scalar multiplication expression
//        (\f$ \vec{a}=(B*s1)*(\vec{c}*s2) \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param vec The right-hand side sparse vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense matrix-scalar multiplication and a sparse vector-scalar multiplication. It
// restructures the expression \f$ \vec{a}=(B*s1)*(\vec{c}*s2) \f$ to the expression
// \f$ \vec{a}=(B*\vec{c})*(s1*s2) \f$.
*/
template< typename MT     // Type of the dense matrix of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool SO         // Storage order of the left-hand side expression
        , typename VT     // Type of the sparse vector of the right-hand side expression
        , typename ST2 >  // Type of the scalar of the right-hand side expression
inline const typename MultExprTrait< DMatScalarMultExpr<MT,ST1,SO>, SVecScalarMultExpr<VT,ST2,false> >::Type
   operator*( const DMatScalarMultExpr<MT,ST1,SO>& mat, const SVecScalarMultExpr<VT,ST2,false>& vec )
{
   return ( mat.leftOperand() * vec.leftOperand() ) * ( mat.rightOperand() * vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse vector-scalar
//        multiplication expression and a dense matrix-scalar multiplication expression
//        (\f$ \vec{a}^T=\vec{b}^T*(C*s1) \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side sparse vector-scalar multiplication.
// \param mat The right-hand side dense matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse vector-scalar multiplication and a dense matrix-scalar multiplication. It
// restructures the expression \f$ \vec{a}=(\vec{b}^T*s1)*(C*s2) \f$ to the expression
// \f$ \vec{a}^T=(\vec{b}^T*C)*(s1*s2) \f$.
*/
template< typename VT   // Type of the sparse vector of the left-hand side expression
        , typename ST1  // Type of the scalar of the left-hand side expression
        , typename MT   // Type of the dense matrix of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO >     // Storage order of the right-hand side expression
inline const typename MultExprTrait< SVecScalarMultExpr<VT,ST1,true>, DMatScalarMultExpr<MT,ST2,SO> >::Type
   operator*( const SVecScalarMultExpr<VT,ST1,true>& vec, const DMatScalarMultExpr<MT,ST2,SO>& mat )
{
   return ( vec.leftOperand() * mat.leftOperand() ) * ( vec.rightOperand() * mat.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar multiplication
//        expression and a dense matrix (\f$ A=(B*s1)*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-scalar multiplication.
// \param rhs The right-hand side dense matrix.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense matrix-scalar multiplication and a dense matrix. It restructures the expression
// \f$ A=(B*s1)*C \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1  // Type of the dense matrix of the left-hand side expression
        , typename ST   // Type of the scalar of the left-hand side expression
        , bool SO1      // Storage order of the left-hand side expression
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline const typename MultExprTrait< DMatScalarMultExpr<MT1,ST,SO1>, MT2 >::Type
   operator*( const DMatScalarMultExpr<MT1,ST,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return ( lhs.leftOperand() * (~rhs) ) * lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix and a dense matrix-
//        scalar multiplication expression (\f$ A=(B*s1)*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix.
// \param rhs The right-hand side dense matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense matrix and a dense matrix-scalar multiplication. It restructures the expression
// \f$ A=B*(C*s1) \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the dense matrix of the right-hand side expression
        , typename ST   // Type of the scalar of the right-hand side expression
        , bool SO2 >    // Storage order of the right-hand side expression
inline const typename MultExprTrait< MT1, DMatScalarMultExpr<MT2,ST,SO2> >::Type
   operator*( const DenseMatrix<MT1,SO1>& lhs, const DMatScalarMultExpr<MT2,ST,SO2>& rhs )
{
   return ( (~lhs) * rhs.leftOperand() ) * rhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of two dense matrix-scalar multiplication
//        expressions (\f$ A=(B*s1)*(C*s2) \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-scalar multiplication.
// \param rhs The right-hand side dense matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of
// two dense matrix-scalar multiplication expressions. It restructures the expression
// \f$ A=(B*s1)*(C*s2) \f$ to the expression \f$ A=(B*C)*(s1*s2) \f$.
*/
template< typename MT1  // Type of the dense matrix of the left-hand side expression
        , typename ST1  // Type of the scalar of the left-hand side expression
        , bool SO1      // Storage order of the left-hand side expression
        , typename MT2  // Type of the right-hand side dense matrix
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO2 >    // Storage order of the right-hand side expression
inline const typename MultExprTrait< DMatScalarMultExpr<MT1,ST1,SO1>, DMatScalarMultExpr<MT2,ST2,SO2> >::Type
   operator*( const DMatScalarMultExpr<MT1,ST1,SO1>& lhs, const DMatScalarMultExpr<MT2,ST2,SO2>& rhs )
{
   return ( lhs.leftOperand() * rhs.leftOperand() ) * ( lhs.rightOperand() * rhs.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar multiplication
//        expression and a sparse matrix (\f$ A=(B*s1)*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-scalar multiplication.
// \param rhs The right-hand side sparse matrix.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of a
// dense matrix-scalar multiplication and a sparse matrix. It restructures the expression
// \f$ A=(B*s1)*C \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1    // Type of the dense matrix of the left-hand side expression
        , typename ST     // Type of the scalar of the left-hand side expression
        , bool SO1        // Storage order of the left-hand side expression
        , typename MT2    // Type of the right-hand side sparse matrix
        , bool SO2 >      // Storage order of the right-hand side sparse matrix
inline const typename MultExprTrait< DMatScalarMultExpr<MT1,ST,SO1>, MT2 >::Type
   operator*( const DMatScalarMultExpr<MT1,ST,SO1>& lhs, const SparseMatrix<MT2,SO2>& rhs )
{
   return ( lhs.leftOperand() * (~rhs) ) * lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse matrix and a dense matrix-
//        scalar multiplication expression (\f$ A=(B*s1)*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix.
// \param rhs The right-hand side dense matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of a
// sparse matrix and a dense matrix-scalar multiplication. It restructures the expression
// \f$ A=B*(C*s1) \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , bool SO1        // Storage order of the left-hand side sparse matrix
        , typename MT2    // Type of the dense matrix of the right-hand side expression
        , typename ST     // Type of the scalar of the right-hand side expression
        , bool SO2 >      // Storage order of the right-hand side expression
inline const typename MultExprTrait< MT1, DMatScalarMultExpr<MT2,ST,SO2> >::Type
   operator*( const SparseMatrix<MT1,SO1>& lhs, const DMatScalarMultExpr<MT2,ST,SO2>& rhs )
{
   return ( (~lhs) * rhs.leftOperand() ) * rhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar
//        multiplication expression and a sparse matrix-scalar multiplication expression
//        (\f$ A=(B*s1)*(C*s2) \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix-scalar multiplication.
// \param vec The right-hand side sparse matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of a dense
// matrix-scalar multiplication and a sparse matrix-scalar multiplication. It restructures the
// expression \f$ A=(B*s1)*(C*s2) \f$ to the expression \f$ A=(B*C)*(s1*s2) \f$.
*/
template< typename MT1  // Type of the dense matrix of the left-hand side expression
        , typename ST1  // Type of the scalar of the left-hand side expression
        , bool SO1      // Storage order of the left-hand side expression
        , typename MT2  // Type of the sparse matrix of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO2 >    // Storage order of the right-hand side expression
inline const typename MultExprTrait< DMatScalarMultExpr<MT1,ST1,SO1>, SMatScalarMultExpr<MT2,ST2,SO2> >::Type
   operator*( const DMatScalarMultExpr<MT1,ST1,SO1>& mat, const SMatScalarMultExpr<MT2,ST2,SO2>& vec )
{
   return ( mat.leftOperand() * vec.leftOperand() ) * ( mat.rightOperand() * vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse matrix-scalar
//        multiplication expression and a dense matrix-scalar multiplication expression
//        (\f$ A=(B*s1)*(C*s2) \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side sparse matrix-scalar multiplication.
// \param vec The right-hand side dense matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of a sparse
// matrix-scalar multiplication and a dense matrix-scalar multiplication. It restructures the
// expression \f$ A=(B*s1)*(C*s2) \f$ to the expression \f$ A=(B*C)*(s1*s2) \f$.
*/
template< typename MT1  // Type of the sparse matrix of the left-hand side expression
        , typename ST1  // Type of the scalar of the left-hand side expression
        , bool SO1      // Storage order of the left-hand side expression
        , typename MT2  // Type of the dense matrix of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO2 >    // Storage order of the right-hand side expression
inline const typename MultExprTrait< SMatScalarMultExpr<MT1,ST1,SO1>, DMatScalarMultExpr<MT2,ST2,SO2> >::Type
   operator*( const SMatScalarMultExpr<MT1,ST1,SO1>& mat, const DMatScalarMultExpr<MT2,ST2,SO2>& vec )
{
   return ( mat.leftOperand() * vec.leftOperand() ) * ( mat.rightOperand() * vec.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATSCALARMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename ST2 >
struct DMatScalarMultTrait< DMatScalarMultExpr<MT,ST1,false>, ST2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<MT,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATSCALARMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename ST2 >
struct TDMatScalarMultTrait< DMatScalarMultExpr<MT,ST1,true>, ST2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<MT,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATSCALARDIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename ST2 >
struct DMatScalarDivTrait< DMatScalarMultExpr<MT,ST1,false>, ST2 >
{
 private:
   //**********************************************************************************************
   enum { condition = IsFloatingPoint<typename MathTrait<ST1,ST2>::DivType>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename DMatScalarMultTrait<MT,typename MathTrait<ST1,ST2>::DivType>::Type  T1;
   typedef typename DMatScalarDivTrait<MT,typename MathTrait<ST1,ST2>::DivType>::Type   T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SelectType<condition,T1,T2>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATSCALARDIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename ST2 >
struct TDMatScalarDivTrait< DMatScalarMultExpr<MT,ST1,true>, ST2 >
{
 private:
   //**********************************************************************************************
   enum { condition = IsFloatingPoint<typename MathTrait<ST1,ST2>::DivType>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   typedef typename TDMatScalarMultTrait<MT,typename MathTrait<ST1,ST2>::DivType>::Type  T1;
   typedef typename TDMatScalarDivTrait<MT,typename MathTrait<ST1,ST2>::DivType>::Type   T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename SelectType<condition,T1,T2>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST, typename VT >
struct DMatDVecMultTrait< DMatScalarMultExpr<MT,ST,false>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value   &&
                                IsDenseVector<VT>::value && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename DMatDVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename VT, typename ST2 >
struct DMatDVecMultTrait< DMatScalarMultExpr<MT,ST1,false>, DVecScalarMultExpr<VT,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value   &&
                                IsDenseVector<VT>::value && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DVecScalarMultTrait<typename DMatDVecMultTrait<MT,VT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATDVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST, typename VT >
struct TDMatDVecMultTrait< DMatScalarMultExpr<MT,ST,true>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsDenseVector<VT>::value && !IsTransposeVector<VT>::value  &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename TDMatDVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename VT, typename ST2 >
struct TDMatDVecMultTrait< DMatScalarMultExpr<MT,ST1,true>, DVecScalarMultExpr<VT,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsDenseVector<VT>::value && !IsTransposeVector<VT>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DVecScalarMultTrait<typename TDMatDVecMultTrait<MT,VT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TDVecDMatMultTrait< VT, DMatScalarMultExpr<MT,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value  &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecDMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename MT, typename ST2 >
struct TDVecDMatMultTrait< DVecScalarMultExpr<VT,ST1,true>, DMatScalarMultExpr<MT,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDVecScalarMultTrait<typename TDVecDMatMultTrait<VT,MT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDVECTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TDVecTDMatMultTrait< VT, DMatScalarMultExpr<MT,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value   &&
                                IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TDVecTDMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename MT, typename ST2 >
struct TDVecTDMatMultTrait< DVecScalarMultExpr<VT,ST1,true>, DMatScalarMultExpr<MT,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value   &&
                                IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDVecScalarMultTrait<typename TDVecTDMatMultTrait<VT,MT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATSVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST, typename VT >
struct DMatSVecMultTrait< DMatScalarMultExpr<MT,ST,false>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value  && IsRowMajorMatrix<MT>::value   &&
                                IsSparseVector<VT>::value && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename DMatSVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename VT, typename ST2 >
struct DMatSVecMultTrait< DMatScalarMultExpr<MT,ST1,false>, SVecScalarMultExpr<VT,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value  && IsRowMajorMatrix<MT>::value   &&
                                IsSparseVector<VT>::value && !IsTransposeVector<VT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DVecScalarMultTrait<typename DMatSVecMultTrait<MT,VT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATSVECMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST, typename VT >
struct TDMatSVecMultTrait< DMatScalarMultExpr<MT,ST,true>, VT >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value  && IsColumnMajorMatrix<MT>::value &&
                                IsSparseVector<VT>::value && !IsTransposeVector<VT>::value  &&
                                IsNumeric<ST>::value
                              , typename DVecScalarMultTrait<typename TDMatSVecMultTrait<MT,VT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename ST1, typename VT, typename ST2 >
struct TDMatSVecMultTrait< DMatScalarMultExpr<MT,ST1,true>, SVecScalarMultExpr<VT,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value  && IsColumnMajorMatrix<MT>::value &&
                                IsSparseVector<VT>::value && !IsTransposeVector<VT>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DVecScalarMultTrait<typename TDMatSVecMultTrait<MT,VT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TSVECDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TSVecDMatMultTrait< VT, DMatScalarMultExpr<MT,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsDenseMatrix<MT>::value  && IsRowMajorMatrix<MT>::value  &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TSVecDMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename MT, typename ST2 >
struct TSVecDMatMultTrait< SVecScalarMultExpr<VT,ST1,true>, DMatScalarMultExpr<MT,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsDenseMatrix<MT>::value  && IsRowMajorMatrix<MT>::value  &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDVecScalarMultTrait<typename TSVecDMatMultTrait<VT,MT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TSVECTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT, typename ST >
struct TSVecTDMatMultTrait< VT, DMatScalarMultExpr<MT,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value && IsTransposeVector<VT>::value   &&
                                IsDenseMatrix<MT>::value  && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST>::value
                              , typename TDVecScalarMultTrait<typename TSVecTDMatMultTrait<VT,MT>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST1, typename MT, typename ST2 >
struct TSVecTDMatMultTrait< SVecScalarMultExpr<VT,ST1,true>, DMatScalarMultExpr<MT,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseVector<VT>::value && IsTransposeVector<VT>::value   &&
                                IsDenseMatrix<MT>::value  && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDVecScalarMultTrait<typename TSVecTDMatMultTrait<VT,MT>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct DMatDMatMultTrait< DMatScalarMultExpr<MT1,ST,false>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DMatDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename ST >
struct DMatDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DMatDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct DMatDMatMultTrait< DMatScalarMultExpr<MT1,ST1,false>, DMatScalarMultExpr<MT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename DMatDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct DMatTDMatMultTrait< DMatScalarMultExpr<MT1,ST,false>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value    &&
                                IsDenseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DMatTDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename ST >
struct DMatTDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value    &&
                                IsDenseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DMatTDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct DMatTDMatMultTrait< DMatScalarMultExpr<MT1,ST1,false>, DMatScalarMultExpr<MT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value    &&
                                IsDenseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename DMatTDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct TDMatDMatMultTrait< DMatScalarMultExpr<MT1,ST,true>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TDMatDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename ST >
struct TDMatDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TDMatDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct TDMatDMatMultTrait< DMatScalarMultExpr<MT1,ST1,true>, DMatScalarMultExpr<MT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<typename TDMatDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct TDMatTDMatMultTrait< DMatScalarMultExpr<MT1,ST,true>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TDMatTDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename ST >
struct TDMatTDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TDMatTDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct TDMatTDMatMultTrait< DMatScalarMultExpr<MT1,ST1,true>, DMatScalarMultExpr<MT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<typename TDMatTDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATSMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct DMatSMatMultTrait< DMatScalarMultExpr<MT1,ST,false>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DMatSMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct DMatSMatMultTrait< DMatScalarMultExpr<MT1,ST1,false>, SMatScalarMultExpr<MT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename DMatSMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DMATTSMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct DMatTSMatMultTrait< DMatScalarMultExpr<MT1,ST,false>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value    &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename DMatTSMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct DMatTSMatMultTrait< DMatScalarMultExpr<MT1,ST1,false>, SMatScalarMultExpr<MT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsRowMajorMatrix<MT1>::value    &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename DMatTSMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATSMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct TDMatSMatMultTrait< DMatScalarMultExpr<MT1,ST,true>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TDMatSMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct TDMatSMatMultTrait< DMatScalarMultExpr<MT1,ST1,true>, SMatScalarMultExpr<MT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<typename TDMatSMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TDMATTSMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct TDMatTSMatMultTrait< DMatScalarMultExpr<MT1,ST,true>, MT2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TDMatTSMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct TDMatTSMatMultTrait< DMatScalarMultExpr<MT1,ST1,true>, SMatScalarMultExpr<MT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT1>::value  && IsColumnMajorMatrix<MT1>::value &&
                                IsSparseMatrix<MT2>::value && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<typename TDMatTSMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SMATDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct SMatDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename SMatDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct SMatDMatMultTrait< SMatScalarMultExpr<MT1,ST1,false>, DMatScalarMultExpr<MT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename SMatDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SMATTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct SMatTDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value    &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename DMatScalarMultTrait<typename SMatTDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct SMatTDMatMultTrait< SMatScalarMultExpr<MT1,ST1,false>, DMatScalarMultExpr<MT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsRowMajorMatrix<MT1>::value    &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<typename SMatTDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TSMATDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct TSMatDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TSMatDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct TSMatDMatMultTrait< SMatScalarMultExpr<MT1,ST1,true>, DMatScalarMultExpr<MT2,ST2,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsRowMajorMatrix<MT2>::value    &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<typename TSMatDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TSMATTDMATMULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST, typename MT2 >
struct TSMatTDMatMultTrait< MT1, DMatScalarMultExpr<MT2,ST,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST>::value
                              , typename TDMatScalarMultTrait<typename TSMatTDMatMultTrait<MT1,MT2>::Type,ST>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename ST1, typename MT2, typename ST2 >
struct TSMatTDMatMultTrait< SMatScalarMultExpr<MT1,ST1,true>, DMatScalarMultExpr<MT2,ST2,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT1>::value && IsColumnMajorMatrix<MT1>::value &&
                                IsDenseMatrix<MT2>::value  && IsColumnMajorMatrix<MT2>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<typename TSMatTDMatMultTrait<MT1,MT2>::Type,typename MathTrait<ST1,ST2>::MultType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
