//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatSMatSubExpr.h
//  \brief Header file for the sparse matrix/sparse matrix subtraction expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATSMATSUBEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATSMATSUBEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/smp/DenseMatrix.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/util/Assert.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATSMATSUBEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-sparse matrix subtractions.
// \ingroup sparse_matrix_expression
//
// The SMatSMatSubExpr class represents the compile time expression for subtractions between
// sparse matrices.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
class SMatSMatSubExpr : public SparseMatrix< SMatSMatSubExpr<MT1,MT2>, false >
                      , private MatMatSubExpr
                      , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT1::ResultType     RT1;  //!< Result type of the left-hand side sparse matrix expression.
   typedef typename MT2::ResultType     RT2;  //!< Result type of the right-hand side sparse matrix expression.
   typedef typename MT1::ReturnType     RN1;  //!< ReturnType type of the left-hand side sparse matrix expression.
   typedef typename MT2::ReturnType     RN2;  //!< ReturnType type of the right-hand side sparse matrix expression.
   typedef typename MT1::CompositeType  CT1;  //!< Composite type of the left-hand side sparse matrix expression.
   typedef typename MT2::CompositeType  CT2;  //!< Composite type of the right-hand side sparse matrix expression.
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
   typedef typename SubExprTrait<RN1,RN2>::Type  ExprReturnType;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SMatSMatSubExpr<MT1,MT2>                    This;           //!< Type of this SMatSMatSubExpr instance.
   typedef typename SubTrait<RT1,RT2>::Type            ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType           OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef const ResultType  CompositeType;

   //! Composite type of the left-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT1>::value, const MT1, const MT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT2>::value, const MT2, const MT2& >::Type  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatSMatSubExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the subtraction expression.
   // \param rhs The right-hand side sparse matrix operand of the subtraction expression.
   */
   explicit inline SMatSMatSubExpr( const MT1& lhs, const MT2& rhs )
      : lhs_( lhs )  // Left-hand side sparse matrix of the subtraction expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the subtraction expression
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
      return lhs_(i,j) - rhs_(i,j);
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

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return lhs_.nonZeros() + rhs_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return lhs_.nonZeros(i) + rhs_.nonZeros(i);
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
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the subtraction expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the subtraction expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-sparse matrix subtraction to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-sparse
   // matrix subtraction expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO>& lhs, const SMatSMatSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      typedef typename RemoveReference<CT2>::Type::ConstIterator  RightIterator;

      assign( ~lhs, rhs.lhs_ );

      if( !IsResizable<typename MT::ElementType>::value ) {
         subAssign( ~lhs, rhs.rhs_ );
      }
      else
      {
         CT2 B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

         BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
         BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
         BLAZE_INTERNAL_ASSERT( B.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
         BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

         for( size_t i=0UL; i<(~lhs).rows(); ++i ) {
            const RightIterator end( B.end(i) );
            for( RightIterator element=B.begin(i); element!=end; ++element ) {
               if( isDefault( (~lhs)(i,element->index()) ) )
                  (~lhs)(i,element->index()) = -element->value();
               else
                  (~lhs)(i,element->index()) -= element->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-sparse matrix subtraction to a row-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-sparse
   // matrix subtraction expression to a row-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,false>& lhs, const SMatSMatSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      typedef typename RemoveReference<CT1>::Type::ConstIterator  LeftIterator;
      typedef typename RemoveReference<CT2>::Type::ConstIterator  RightIterator;

      CT1 A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).columns()  , "Invalid number of columns" );

      // Final memory allocation (based on the evaluated operands)
      (~lhs).reserve( A.nonZeros() + B.nonZeros() );

      // Performing the matrix subtraction
      for( size_t i=0UL; i<(~lhs).rows(); ++i )
      {
         const LeftIterator  lend( A.end(i) );
         const RightIterator rend( B.end(i) );

         LeftIterator  l( A.begin(i) );
         RightIterator r( B.begin(i) );

         while( l != lend && r != rend )
         {
            if( l->index() < r->index() ) {
               (~lhs).append( i, l->index(), l->value() );
               ++l;
            }
            else if( l->index() > r->index() ) {
               (~lhs).append( i, r->index(), -r->value() );
               ++r;
            }
            else {
               (~lhs).append( i, l->index(), l->value()-r->value() );
               ++l;
               ++r;
            }
         }

         while( l != lend ) {
            (~lhs).append( i, l->index(), l->value() );
            ++l;
         }

         while( r != rend ) {
            (~lhs).append( i, r->index(), -r->value() );
            ++r;
         }

         (~lhs).finalize( i );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-sparse matrix subtraction to a column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-sparse
   // matrix subtraction expression to a column-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,true>& lhs, const SMatSMatSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      typedef typename RemoveReference<CT1>::Type::ConstIterator  LeftIterator;
      typedef typename RemoveReference<CT2>::Type::ConstIterator  RightIterator;

      CT1 A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).columns()  , "Invalid number of columns" );

      const size_t m( rhs.rows()    );
      const size_t n( rhs.columns() );

      // Counting the number of elements per column
      std::vector<size_t> nonzeros( n, 0UL );
      for( size_t i=0UL; i<m; ++i )
      {
         const LeftIterator  lend( A.end(i) );
         const RightIterator rend( B.end(i) );

         LeftIterator  l( A.begin(i) );
         RightIterator r( B.begin(i) );

         while( l != lend && r != rend )
         {
            if( l->index() < r->index() ) {
               ++nonzeros[l->index()];
               ++l;
            }
            else if( l->index() > r->index() ) {
               ++nonzeros[r->index()];
               ++r;
            }
            else {
               ++nonzeros[l->index()];
               ++l;
               ++r;
            }
         }

         while( l != lend ) {
            ++nonzeros[l->index()];
            ++l;
         }

         while( r != rend ) {
            ++nonzeros[r->index()];
            ++r;
         }
      }

      // Resizing the left-hand side sparse matrix
      for( size_t j=0UL; j<n; ++j ) {
         (~lhs).reserve( j, nonzeros[j] );
      }

      // Performing the matrix subtraction
      for( size_t i=0UL; i<m; ++i )
      {
         const LeftIterator  lend( A.end(i) );
         const RightIterator rend( B.end(i) );

         LeftIterator  l( A.begin(i) );
         RightIterator r( B.begin(i) );

         while( l != lend && r != rend )
         {
            if( l->index() < r->index() ) {
               (~lhs).append( i, l->index(), l->value() );
               ++l;
            }
            else if( l->index() > r->index() ) {
               (~lhs).append( i, r->index(), -r->value() );
               ++r;
            }
            else {
               (~lhs).append( i, l->index(), l->value()-r->value() );
               ++l;
               ++r;
            }
         }

         while( l != lend ) {
            (~lhs).append( i, l->index(), l->value() );
            ++l;
         }

         while( r != rend ) {
            (~lhs).append( i, r->index(), -r->value() );
            ++r;
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-sparse matrix subtraction to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side subtraction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // sparse matrix subtraction expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO>& lhs, const SMatSMatSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAddAssign( ~lhs, rhs.lhs_ );
      smpSubAssign( ~lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-sparse matrix subtraction to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side subtraction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // sparse matrix subtraction expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO>& lhs, const SMatSMatSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSubAssign( ~lhs, rhs.lhs_ );
      smpAddAssign( ~lhs, rhs.rhs_ );
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
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
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
/*!\brief Subtraction operator for the subtraction of two row-major sparse matrices (\f$ A=B-C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the matrix subtraction.
// \param rhs The right-hand side sparse matrix to be subtracted from the left-hand side matrix.
// \return The difference of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match
//
// This operator represents the subtraction of two row-major sparse matrices:

   \code
   using blaze::rowMajor;

   blaze::CompressedMatrix<double,rowMajor> A, B, C;
   // ... Resizing and initialization
   C = A - B;
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the two involved matrix element types \a T1::ElementType and \a T2::ElementType.
// Both matrix types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the SubTrait class template.\n
// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side sparse matrix
inline const SMatSMatSubExpr<T1,T2>
   operator-( const SparseMatrix<T1,false>& lhs, const SparseMatrix<T2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   return SMatSMatSubExpr<T1,T2>( ~lhs, ~rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool AF >
struct SubmatrixExprTrait< SMatSMatSubExpr<MT1,MT2>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename SubExprTrait< typename SubmatrixExprTrait<const MT1,AF>::Type
                                , typename SubmatrixExprTrait<const MT2,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct RowExprTrait< SMatSMatSubExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SubExprTrait< typename RowExprTrait<const MT1>::Type
                                , typename RowExprTrait<const MT2>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct ColumnExprTrait< SMatSMatSubExpr<MT1,MT2> >
{
 public:
   //**********************************************************************************************
   typedef typename SubExprTrait< typename ColumnExprTrait<const MT1>::Type
                                , typename ColumnExprTrait<const MT2>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
