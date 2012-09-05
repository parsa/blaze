//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatScalarDivExpr.h
//  \brief Header file for the dense matrix/scalar division expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSCALARDIVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSCALARDIVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Expression.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATSCALARDIVEXPRHELPER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Helper class for divisions of a dense matrix by a scalar.
// \ingroup dense_matrix_expression
//
// The DMatScalarDivExprHelper class is an auxiliary class to define the return type of the
// division between a dense matrix and a scalar value.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar value
        , bool SO >    // Storage order
struct DMatScalarDivExprHelper
{
 public:
   //**Type definitions****************************************************************************
   //! Scalar type for the instantiation of the resulting expression object.
   typedef typename MathTrait< typename BaseElementType<MT>::Type, ST >::DivType  ScalarType;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation of the dense matrix/scalar division return type.
   enum { value = IsFloatingPoint<ScalarType>::value };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Resulting type of the division between the given dense matrix and scalar value.
   typedef typename SelectType< value,
                                DMatScalarMultExpr<MT,ScalarType,SO>,
                                DMatScalarDivExpr<MT,ScalarType,SO> >::Type  Type;
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ScalarType );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DMATSCALARDIVEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for divisions of a dense matrix by a scalar.
// \ingroup dense_matrix_expression
//
// The DMatScalarDivExpr class represents the compile time expression for divisions of dense
// matrices and by scalar values.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar value
        , bool SO >    // Storage order
class DMatScalarDivExpr : public DenseMatrix< DMatScalarDivExpr<MT,ST,SO>, SO >
                        , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     RT;  //!< Result type of the dense matrix expression.
   typedef typename MT::ReturnType     RN;  //!< Return type of the dense matrix expression.
   typedef typename MT::CompositeType  CT;  //!< Composite type of the dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the division expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the division expression. In case the dense matrix operand
       requires an intermediate evaluation, \a useAssign will be set to \a true and the
       division expression will be evaluated via the \a assign function family. Otherwise
       \a useAssign will be set to \a false and the expression will be evaluated via the
       subscript operator. */
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
   typedef DMatScalarDivExpr<MT,ST,SO>               This;           //!< Type of this DMatScalarDivExpr instance.
   typedef typename MathTrait<RT,ST>::DivType        ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType         OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType        TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType          ElementType;    //!< Resulting element type.
   typedef const typename DivExprTrait<RN,ST>::Type  ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DMatScalarDivExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  LeftOperand;

   //! Composite type of the right-hand side scalar value.
   typedef typename MathTrait< typename BaseElementType<MT>::Type, ST >::DivType  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<MT>::value };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatScalarDivExpr class.
   //
   // \param matrix The left-hand side dense matrix of the division expression.
   // \param scalar The right-hand side scalar of the division expression.
   */
   explicit inline DMatScalarDivExpr( const MT& matrix, ST scalar )
      : matrix_( matrix )  // Left-hand side dense matrix of the division expression
      , scalar_( scalar )  // Right-hand side scalar of the division expression
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
      return matrix_(i,j) / scalar_;
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
   LeftOperand  matrix_;  //!< Left-hand side dense matrix of the division expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the division expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-scalar division to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // division expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( DenseMatrix<MT2,false>& lhs, const DMatScalarDivExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );
      (~lhs) /= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-scalar division to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // division expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( DenseMatrix<MT2,true>& lhs, const DMatScalarDivExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );
      (~lhs) /= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-scalar division to a row-major sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // division expression to a sparse matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target sparse matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( SparseMatrix<MT2,false>& lhs, const DMatScalarDivExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );

      for( size_t i=0UL; i<(~lhs).rows(); ++i )
      {
         typename MT2::Iterator element( (~lhs).begin(i) );
         const typename MT2::Iterator end( (~lhs).end(i) );

         for( ; element!=end; ++element )
            element->value() /= rhs.scalar_;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-scalar division to a column-major sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-scalar
   // division expression to a sparse matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2 >  // Type of the target sparse matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( SparseMatrix<MT2,true>& lhs, const DMatScalarDivExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.matrix_ );

      for( size_t j=0UL; j<(~lhs).columns(); ++j )
      {
         typename MT2::Iterator element( (~lhs).begin(j) );
         const typename MT2::Iterator end( (~lhs).end(j) );

         for( ; element!=end; ++element )
            element->value() /= rhs.scalar_;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-scalar division to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // scalar division expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the matrix
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatScalarDivExpr& rhs )
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
   /*!\brief Subtraction assignment of a dense matrix-scalar division to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix-
   // scalar division expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the matrix operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatScalarDivExpr& rhs )
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
   BLAZE_CONSTRAINT_MUST_NOT_BE_FLOATING_POINT_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_NOT_BE_FLOATING_POINT_TYPE( ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
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
/*!\brief Division operator for the division of a dense matrix by a scalar value (\f$ A=B/s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This operator represents the division of a dense matrix by a scalar value:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = A / 0.24;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a T1::ElementType and \a T2. Note that this operator only
// works for scalar values of built-in data type.
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename T1    // Type of the left-hand side dense matrix
        , bool SO        // Storage order of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsNumeric<T2>
                              , typename DMatScalarDivExprHelper<T1,T2,SO>::Type >::Type
   operator/( const DenseMatrix<T1,SO>& mat, T2 scalar )
{
   BLAZE_USER_ASSERT( scalar != T2(0), "Division by zero detected" );

   typedef DMatScalarDivExprHelper<T1,T2,SO>  Helper;
   typedef typename Helper::ScalarType        ScalarType;

   if( Helper::value ) {
      return typename Helper::Type( ~mat, ScalarType(1)/ScalarType(scalar) );
   }
   else {
      return typename Helper::Type( ~mat, scalar );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix-scalar division
//        expression and a scalar value (\f$ A=(B/s1)*s2 \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix-scalar division.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the multiplication of a
// dense matrix-scalar division expression and a scalar value.
*/
template< typename MT     // Type of the dense matrix of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool SO         // Storage order of the dense matrix
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsFloatingPoint< typename MathTrait<ST2,ST1>::DivType >
                              , typename MultExprTrait< DMatScalarDivExpr<MT,ST1,SO>, ST2 >::Type >::Type
   operator*( const DMatScalarDivExpr<MT,ST1,SO>& mat, ST2 scalar )
{
   return mat.leftOperand() * ( scalar / mat.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense matrix-
//        scalar division expression (\f$ A=s2*(B/s1) \f$).
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param mat The right-hand side dense matrix-scalar division.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a dense matrix-scalar division expression.
*/
template< typename ST1  // Type of the left-hand side scalar
        , typename MT   // Type of the dense matrix of the right-hand side expression
        , typename ST2  // Type of the scalar of the right-hand side expression
        , bool SO >     // Storage order of the dense matrix
inline const typename EnableIf< IsFloatingPoint< typename MathTrait<ST1,ST2>::DivType >
                              , typename MultExprTrait< ST1, DMatScalarDivExpr<MT,ST2,SO> >::Type >::Type
   operator*( ST1 scalar, const DMatScalarDivExpr<MT,ST2,SO>& mat )
{
   return mat.leftOperand() * ( scalar / mat.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a dense matrix-scalar division expression
//        and a scalar value (\f$ A=(B/s1)/s2 \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix-scalar division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the division of a dense
// matrix-scalar division expression and a scalar value.
*/
template< typename MT     // Type of the dense matrix of the left-hand side expression
        , typename ST1    // Type of the scalar of the left-hand side expression
        , bool SO         // Storage order of the dense matrix
        , typename ST2 >  // Type of the right-hand side scalar
inline const typename EnableIf< IsNumeric<ST2>
                              , typename DMatScalarDivExprHelper<MT,typename MathTrait<ST1,ST2>::MultType,SO>::Type >::Type
   operator/( const DMatScalarDivExpr<MT,ST1,SO>& mat, ST2 scalar )
{
   BLAZE_USER_ASSERT( scalar != ST2(0), "Division by zero detected" );

   typedef typename MathTrait<ST1,ST2>::MultType    MultType;
   typedef DMatScalarDivExprHelper<MT,MultType,SO>  Helper;

   if( Helper::value ) {
      return typename Helper::Type( mat.leftOperand(), MultType(1)/( mat.rightOperand() * scalar ) );
   }
   else {
      return typename Helper::Type( mat.leftOperand(), mat.rightOperand() * scalar );
   }
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
struct DMatScalarMultTrait< DMatScalarDivExpr<MT,ST1,false>, ST2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename DMatScalarMultTrait<MT,typename MathTrait<ST2,ST1>::DivType>::Type
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
struct TDMatScalarMultTrait< DMatScalarDivExpr<MT,ST1,true>, ST2 >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsDenseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value &&
                                IsNumeric<ST1>::value && IsNumeric<ST2>::value
                              , typename TDMatScalarMultTrait<MT,typename MathTrait<ST2,ST1>::DivType>::Type
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
