//=================================================================================================
/*!
//  \file blaze/math/expressions/TDVecTSMatMultExpr.h
//  \brief Header file for the transpose dense vector/transpose sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDVECTSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDVECTSMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <boost/type_traits/remove_reference.hpp>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Expression.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDVECSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense vector-transpose sparse matrix multiplications.
// \ingroup dense_vector_expression
//
// The TDVecTSMatMultExpr class represents the compile time expression for multiplications
// between transpose dense vectors and column-major sparse matrices.
*/
template< typename VT    // Type of the left-hand side dense vector
        , typename MT >  // Type of the right-hand side sparse matrix
class TDVecTSMatMultExpr : public DenseVector< TDVecTSMatMultExpr<VT,MT>, true >
                         , private Expression
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     VRT;  //!< Result type of the left-hand side dense vector expression.
   typedef typename MT::ResultType     MRT;  //!< Result type of the right-hand side sparse matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the left-hand side dense vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case the sparse matrix
       expression requires an intermediate evaluation or the dense vector expression is a
       compound expression, \a useAssign will be set to \a true and the addition expression
       will be evaluated via the \a assign function family. Otherwise \a useAssign will be
       set to \a false and the expression will be evaluated via the subscript operator. */
   enum { useAssign = ( IsExpression<VT>::value || !IsReference<MCT>::value ) };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   struct UseAssign {
      enum { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TDVecTSMatMultExpr<VT,MT>              This;           //!< Type of this TDVecTSMatMultExpr instance.
   typedef typename MathTrait<VRT,MRT>::MultType  ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType     TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType       ElementType;    //!< Resulting element type.
   typedef const ElementType                      ReturnType;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const TDVecTSMatMultExpr& >::Type  CompositeType;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  LeftOperand;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  RightOperand;

   //! Composite type of the left-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VRT, VCT >::Type  LT;

   //! Composite type of the right-hand side sparse matrix expression.
   typedef MCT  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = ( !IsExpression<VT>::value ) ||
                     ( IsReference<MCT>::value && ( !IsExpression<MT>::value || CanAlias<MT>::value ) ) };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDVecTSMatMultExpr class.
   */
   explicit inline TDVecTSMatMultExpr( const VT& vec, const MT& mat )
      : vec_( vec )  // Left-hand side dense vector of the multiplication expression
      , mat_( mat )  // Right-hand side sparse matrix of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( vec_.size() == mat.rows(), "Invalid vector and matrix sizes" );
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

      typedef typename boost::remove_reference<MCT>::type::ConstIterator  ConstIterator;

      VCT x( vec_ );  // Evaluation of the left-hand side sparse vector operand
      MCT A( mat_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( x.size()    == vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == mat_.columns(), "Invalid number of columns" );

      const ConstIterator end( A.end(index) );
      ConstIterator element( A.begin(index) );
      ElementType res;

      if( element != end ) {
         res = x[element->index()] * element->value();
         ++element;
         for( ; element!=end; ++element )
            res += x[element->index()] * element->value();
      }
      else {
         reset( res );
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

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
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
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the given alias is contained in this expression, \a false if not.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return ( !IsExpression<VT>::value && vec_.isAliased( alias ) ) ||
             ( IsReference<MCT>::value && mat_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vec_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand mat_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a transpose dense vector-transpose sparse matrix multiplication to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense vector-
   // transpose sparse matrix multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler in
   // case either the left-hand side vector operand is a compound expression or the right-hand
   // side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<RT>::type::ConstIterator  ConstIterator;

      if( rhs.mat_.rows() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT x( rhs.vec_ );  // Evaluation of the left-hand side dense vector operator
      RT A( rhs.mat_ );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         const ConstIterator end( A.end(j) );
         ConstIterator element( A.begin(j) );

         if( element == end ) {
            reset( (~lhs)[j] );
            continue;
         }

         (~lhs)[j] = x[element->index()] * element->value();
         ++element;
         for( ; element!=end; ++element ) {
            (~lhs)[j] += x[element->index()] * element->value();
         }
      }
   }
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a transpose dense vector-transpose sparse matrix multiplication to a
   //        sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense vector-
   // transpose sparse matrix multiplication expression to a sparse vector. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler in
   // case either the left-hand side vector operand is a compound expression or the right-hand
   // side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( SparseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      assign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a transpose dense vector-transpose sparse matrix multiplication
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due
   // to the explicit application of the SFINAE principle, this operator can only be selected
   // by the compiler in case either the left-hand side vector operand is a compound expression
   // or the right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      addAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<RT>::type::ConstIterator  ConstIterator;

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }

      LT x( rhs.vec_ );  // Evaluation of the left-hand side dense vector operator
      RT A( rhs.mat_ );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j ) {
         const ConstIterator end( A.end(j) );
         ConstIterator element( A.begin(j) );

         for( ; element!=end; ++element ) {
            (~lhs)[j] += x[element->index()] * element->value();
         }
      }
   }
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this operator can only be selected by
   // the compiler in case either the left-hand side vector operand is a compound expression
   // or the right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      subAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<RT>::type::ConstIterator  ConstIterator;

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }

      LT x( rhs.vec_ );  // Evaluation of the left-hand side dense vector operator
      RT A( rhs.mat_ );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j ) {
         const ConstIterator end( A.end(j) );
         ConstIterator element( A.begin(j) );

         for( ; element!=end; ++element ) {
            (~lhs)[j] -= x[element->index()] * element->value();
         }
      }
   }
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this operator can only be selected by the
   // compiler in case either the left-hand side vector operand is a compound expression or the
   // right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      multAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      multAssign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE( VT );
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
/*!\brief Multiplication operator for the multiplication of a transpose dense vector and a
//        column-major sparse matrix (\f$ \vec{y}^T=\vec{x}^T*A \f$).
// \ingroup sparse_matrix
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side column-major sparse matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose dense vector and a column-major
// sparse matrix:

   \code
   using blaze::rowVector;
   using blaze::columnMajor;

   blaze::DynamicVector<double,rowVector> x, y;
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// The operator returns an expression representing a transpose dense vector of the higher-order
// element type of the two involved element types \a T1::ElementType and \a T2::ElementType.
// Both the dense matrix type \a T1 and the dense vector type \a T2 as well as the two element
// types \a T1::ElementType and \a T2::ElementType have to be supported by the MathTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of rows of
// the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side dense vector
        , typename T2 >  // Type of the right-hand side sparse matrix
inline const typename DisableIf< IsMatMatMultExpr<T2>, TDVecTSMatMultExpr<T1,T2> >::Type
   operator*( const DenseVector<T1,true>& vec, const SparseMatrix<T2,true>& mat )
{
   if( (~vec).size() != (~mat).rows() )
      throw std::invalid_argument( "Vector and matrix sizes do not match" );

   return TDVecTSMatMultExpr<T1,T2>( ~vec, ~mat );
}
//*************************************************************************************************

} // namespace blaze

#endif
