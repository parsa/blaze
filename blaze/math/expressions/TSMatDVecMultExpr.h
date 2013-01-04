//=================================================================================================
/*!
//  \file blaze/math/expressions/TSMatDVecMultExpr.h
//  \brief Header file for the transpose sparse matrix/dense vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TSMATDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSMATDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <boost/type_traits/remove_reference.hpp>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TSMATDVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose sparse matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// The TSMatDVecMultExpr class represents the compile time expression for multiplications
// between column-major sparse matrices and dense vectors.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , typename VT >  // Type of the right-hand side dense vector
class TSMatDVecMultExpr : public DenseVector< TSMatDVecMultExpr<MT,VT>, false >
                        , private Expression
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     MRT;  //!< Result type of the left-hand side sparse matrix expression.
   typedef typename VT::ResultType     VRT;  //!< Result type of the right-hand side dense vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the left-hand side sparse matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TSMatDVecMultExpr<MT,VT>            This;           //!< Type of this TSMatDVecMultExpr instance.
   typedef typename MultTrait<MRT,VRT>::Type   ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.
   typedef const ElementType                   ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                    CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side sparse matrix operand.
   typedef MCT  LT;

   //! Type for the assignment of the right-hand side dense matrix operand.
   typedef typename SelectType< IsComputation<VT>::value, const VRT, VCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = ( IsComputation<MT>::value && !RequiresEvaluation<MT>::value && CanAlias<MT>::value ) ||
                     ( !IsComputation<VT>::value ) };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TSMatDVecMultExpr class.
   //
   // \param mat The left-hand side sparse matrix operand of the multiplication expression.
   // \param vec The right-hand side dense vector operand of the multiplication expression.
   */
   explicit inline TSMatDVecMultExpr( const MT& mat, const VT& vec )
      : mat_( mat )  // Left-hand side sparse matrix of the multiplication expression
      , vec_( vec )  // Right-hand side dense vector of the multiplication expression
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

      ElementType res;

      if( vec_.size() > 0UL ) {
         res = mat_(index,0UL) * vec_[0UL];
         for( size_t j=1UL; j<vec_.size(); ++j )
            res += mat_(index,j) * vec_[j];
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
      return mat_.rows();
   }
   //**********************************************************************************************

   //**Left function*******************************************************************************
   /*!\brief Returns the left-hand side transpose sparse matrix operand.
   //
   // \return The left-hand side transpose sparse matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return mat_;
   }
   //**********************************************************************************************

   //**Right function******************************************************************************
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const {
      return vec_;
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
      return ( IsComputation<MT>::value && !RequiresEvaluation<MT>::value &&
               CanAlias<MT>::value && mat_.isAliased( alias ) ) ||
             ( !IsComputation<VT>::value && vec_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side sparse matrix of the multiplication expression.
   RightOperand vec_;  //!< Right-hand side dense vector of the multiplication expression.
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose sparse matrix-dense vector multiplication to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the default assignment of a transpose sparse matrix-dense vector
   // multiplication expression to a dense vector. This assign function is used in case the
   // element type of the target vector is resizable.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename EnableIf< IsResizable<typename VT1::ElementType> >::Type
      assign( DenseVector<VT1,false>& lhs, const TSMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<LT>::type::ConstIterator  ConstIterator;

      reset( ~lhs );

      if( rhs.mat_.columns() == 0UL ) return;

      LT A( rhs.mat_ );  // Evaluation of the left-hand side sparse matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         ConstIterator element( A.begin(j) );
         const ConstIterator end( A.end(j) );

         for( ; element!=end; ++element ) {
            if( isDefault( (~lhs)[element->index()] ) )
               (~lhs)[element->index()] = element->value() * x[j];
            else
               (~lhs)[element->index()] += element->value() * x[j];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a transpose sparse matrix-dense vector multiplication to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // dense vector multiplication expression to a dense vector. This assign function is used in
   // case the element type of the target vector is not resizable.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline typename DisableIf< IsResizable<typename VT1::ElementType> >::Type
      assign( DenseVector<VT1,false>& lhs, const TSMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<LT>::type::ConstIterator  ConstIterator;

      reset( ~lhs );

      if( rhs.mat_.columns() == 0UL ) return;

      LT A( rhs.mat_ );  // Evaluation of the left-hand side sparse matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         ConstIterator element( A.begin(j) );
         const ConstIterator end( A.end(j) );

         for( ; element!=end; ++element ) {
            (~lhs)[element->index()] += element->value() * x[j];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-dense vector multiplication to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const TSMatDVecMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-dense vector multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // dense vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this operator can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const TSMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<LT>::type::ConstIterator  ConstIterator;

      if( rhs.mat_.columns() == 0UL ) return;

      LT A( rhs.mat_ );  // Evaluation of the left-hand side sparse matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         ConstIterator element( A.begin(j) );
         const ConstIterator end( A.end(j) );

         for( ; element!=end; ++element ) {
            (~lhs)[element->index()] += element->value() * x[j];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-dense vector multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // dense vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this operator can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const TSMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename boost::remove_reference<LT>::type::ConstIterator  ConstIterator;

      if( rhs.mat_.columns() == 0UL ) return;

      LT A( rhs.mat_ );  // Evaluation of the left-hand side sparse matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         ConstIterator element( A.begin(j) );
         const ConstIterator end( A.end(j) );

         for( ; element!=end; ++element ) {
            (~lhs)[element->index()] -= element->value() * x[j];
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse matrix-dense vector multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this operator can only be selected by the compiler
   // in case either the left-hand side matrix operand requires an intermediate evaluation or
   // the right-hand side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const TSMatDVecMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      multAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT );
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
/*!\brief Multiplication operator for the multiplication of a column-major sparse matrix and a
//        dense vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side column-major sparse matrix for the multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a column-major sparse matrix and a dense
// vector:

   \code
   using blaze::columnMajor;
   using blaze::columnVector;

   blaze::CompressedMatrix<double,columnMajor> A;
   blaze::DynamicVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved element types \a T1::ElementType and \a T2::ElementType. Both the
// sparse matrix type \a T1 and the dense vector type \a T2 as well as the two element types
// \a T1::ElementType and \a T2::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side dense vector
inline const typename DisableIf< IsMatMatMultExpr<T1>, TSMatDVecMultExpr<T1,T2> >::Type
   operator*( const SparseMatrix<T1,true>& mat, const DenseVector<T2,false>& vec )
{
   if( (~mat).columns() != (~vec).size() )
      throw std::invalid_argument( "Matrix and vector sizes do not match" );

   return TSMatDVecMultExpr<T1,T2>( ~mat, ~vec );
}
//*************************************************************************************************

} // namespace blaze

#endif
