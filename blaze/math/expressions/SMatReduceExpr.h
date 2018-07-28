//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatReduceExpr.h
//  \brief Header file for the sparse matrix reduce expression
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATREDUCEEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATREDUCEEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatReduceExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/functors/Add.h>
#include <blaze/math/functors/Mult.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/HasMember.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for row-major sparse matrix reduction operations.
// \ingroup dense_vector_expression
//
// The SMatReduceExpr class represents the compile time expression for partial reduction operations
// of row-major sparse matrices.
*/
template< size_t RF      // Reduction flag
        , typename MT    // Type of the sparse matrix
        , typename OP >  // Type of the reduction operation
struct SMatReduceExpr
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-WISE REDUCTION OPERATIONS OF ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for column-wise row-major sparse matrix reduction operations.
// \ingroup dense_vector_expression
//
// This specialization of the SMatReduceExpr class template represents the compile time expression
// for column-wise reduction operations of row-major sparse matrices.
*/
template< typename MT    // Type of the sparse matrix
        , typename OP >  // Type of the reduction operation
struct SMatReduceExpr<0UL,MT,OP>
   : public MatReduceExpr<0UL, DenseVector< SMatReduceExpr<0UL,MT,OP>, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;     //!< Result type of the sparse matrix expression.
   using OT = OppositeType_t<MT>;   //!< Opposite type of the sparse matrix expression.
   using ET = ElementType_t<MT>;    //!< Element type of the sparse matrix expression.
   using CT = CompositeType_t<MT>;  //!< Composite type of the sparse matrix expression.
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the sparse matrix operand is not SMP assignable and requires an intermediate
       evaluation, the variable is set to 1 and the expression specific evaluation strategy is
       selected. Otherwise the variable is set to 0 and the default strategy is chosen. */
   template< typename VT >
   static constexpr bool UseSMPAssign_v = ( !MT::smpAssignable && RequiresEvaluation_v<MT> );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This          = SMatReduceExpr<0UL,MT,OP>;    //!< Type of this SMatReduceExpr instance.
   using ResultType    = ReduceTrait_t<RT,OP,0UL>;     //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< Resulting SIMD element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Data type of the custom unary operation.
   using Operation = OP;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatReduceExpr class.
   //
   // \param sm The matrix operand of the reduction expression.
   */
   explicit inline SMatReduceExpr( const MT& sm, OP op ) noexcept
      : sm_( sm )  // Sparse matrix of the reduction expression
      , op_( op )  // The reduction operation
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < sm_.columns(), "Invalid vector access index" );
      return reduce( column( sm_, index, unchecked ), op_ );
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= sm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return sm_.columns();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the sparse matrix operand.
   //
   // \return The sparse matrix operand.
   */
   inline Operand operand() const noexcept {
      return sm_;
   }
   //**********************************************************************************************

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the reduction operation.
   //
   // \return A copy of the reduction operation.
   */
   inline Operation operation() const {
      return op_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an aliasing effect is possible, \a false if not.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( sm_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the given alias is contained in this expression, \a false if not.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return ( sm_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return false;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return sm_.canSMPAssign() || ( size() > SMP_SMATDVECMULT_THRESHOLD );  // TODO: rethink
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand   sm_;  //!< Sparse matrix of the reduction expression.
   Operation op_;  //!< The reduction operation.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a column-wise row-major sparse matrix reduction operation to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a column-wise row-major
   // sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      assign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a column-wise row-major sparse matrix reduction operation to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a column-wise row-major
   // sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( SparseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a column-wise row-major sparse matrix reduction operation to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      addAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( SparseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      addAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      subAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( SparseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      subAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      multAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( SparseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      divAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a column-wise row-major sparse matrix reduction operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side reduction expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a column-wise
   // row-major sparse matrix reduction expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void divAssign( SparseVector<VT1,true>& lhs, const SMatReduceExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      divAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to vectors*******************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a column-wise row-major sparse matrix reduction operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a column-wise row-major
   // sparse matrix reduction expression to a vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression specific
   // parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpAssign( Vector<VT1,true>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to vectors**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a column-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a column-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpAddAssign( Vector<VT1,true>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpAddAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a column-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a column-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpSubAssign( Vector<VT1,true>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpSubAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a column-wise row-major sparse matrix reduction
   //        operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // column-wise row-major sparse matrix reduction expression to a vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpMultAssign( Vector<VT1,true>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpMultAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to vectors**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a column-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a column-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpDivAssign( Vector<VT1,true>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpDivAssign( ~lhs, reduce<0UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ROW-WISE REDUCTION OPERATIONS OF ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for row-wise row-major sparse matrix reduction operations.
// \ingroup dense_vector_expression
//
// This specialization of the SMatReduceExpr class template represents the compile time expression
// for row-wise reduction operations of row-major sparse matrices.
*/
template< typename MT    // Type of the sparse matrix
        , typename OP >  // Type of the reduction operation
struct SMatReduceExpr<1UL,MT,OP>
   : public MatReduceExpr<1UL, DenseVector< SMatReduceExpr<1UL,MT,OP>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;   //!< Result type of the sparse matrix expression.
   using ET = ElementType_t<MT>;  //!< Element type of the sparse matrix expression.
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the reduction expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the reduction expression. In case the sparse matrix
       operand requires an intermediate evaluation, \a useAssign will be set to 1 and the
       reduction expression will be evaluated via the \a assign function family. Otherwise
       \a useAssign will be set to 0 and the expression will be evaluated via the subscript
       operator. */
   static constexpr bool useAssign = RequiresEvaluation_v<MT>;

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the sparse matrix operand is not SMP assignable and requires an intermediate
       evaluation, the variable is set to 1 and the expression specific evaluation strategy is
       selected. Otherwise the variable is set to 0 and the default strategy is chosen. */
   template< typename VT >
   static constexpr bool UseSMPAssign_v = ( !MT::smpAssignable && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This          = SMatReduceExpr<1UL,MT,OP>;    //!< Type of this SMatReduceExpr instance.
   using ResultType    = ReduceTrait_t<RT,OP,1UL>;     //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< Resulting SIMD element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const SMatReduceExpr& >;

   //! Composite type of the left-hand side sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Data type of the custom unary operation.
   using Operation = OP;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::random_access_iterator_tag;  //!< The iterator category.
      using ValueType        = ElementType;                      //!< Type of the underlying elements.
      using PointerType      = ElementType*;                     //!< Pointer return type.
      using ReferenceType    = ElementType&;                     //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                        //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param sm The sparse matrix operand of the reduction expression.
      // \param index Index to the initial matrix row.
      // \param op The reduction operation.
      */
      explicit inline ConstIterator( Operand sm, size_t index, OP op )
         : sm_   ( sm    )  // Sparse matrix of the reduction expression
         , index_( index )  // Index to the current matrix row
         , op_   ( op    )  // The reduction operation
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ConstIterator& operator+=( size_t inc ) {
         index_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline ConstIterator& operator-=( size_t dec ) {
         index_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ConstIterator& operator++() {
         ++index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator++( int ) {
         return ConstIterator( index_++ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ConstIterator& operator--() {
         --index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator--( int ) {
         return ConstIterator( index_-- );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReturnType operator*() const {
         return reduce( row( sm_, index_, unchecked ), op_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ConstIterator& rhs ) const {
         return index_ < rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ConstIterator& rhs ) const {
         return index_ > rhs.index_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ConstIterator& rhs ) const {
         return index_ <= rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ConstIterator& rhs ) const {
         return index_ >= rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return index_ - rhs.index_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ConstIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( const ConstIterator& it, size_t inc ) {
         return ConstIterator( it.index_ + inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ConstIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( size_t inc, const ConstIterator& it ) {
         return ConstIterator( it.index_ + inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ConstIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ConstIterator operator-( const ConstIterator& it, size_t dec ) {
         return ConstIterator( it.index_ - dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      Operand sm_;     //!< Sparse matrix of the reduction expression.
      size_t  index_;  //!< Index to the current matrix row.
      OP      op_;     //!< The reduction operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatReduceExpr class.
   //
   // \param sm The matrix operand of the reduction expression.
   // \param op the reduction operation.
   */
   explicit inline SMatReduceExpr( const MT& sm, OP op ) noexcept
      : sm_( sm )  // Sparse matrix of the reduction expression
      , op_( op )  // The reduction operation
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < sm_.rows(), "Invalid vector access index" );
      return reduce( row( sm_, index, unchecked ), op_ );
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= sm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( sm_, 0UL, op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( sm_, size(), op_ );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return sm_.rows();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the sparse matrix operand.
   //
   // \return The sparse matrix operand.
   */
   inline Operand operand() const noexcept {
      return sm_;
   }
   //**********************************************************************************************

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the reduction operation.
   //
   // \return A copy of the reduction operation.
   */
   inline Operation operation() const {
      return op_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an aliasing effect is possible, \a false if not.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( sm_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the given alias is contained in this expression, \a false if not.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return ( sm_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return false;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return sm_.canSMPAssign() || ( size() > SMP_SMATDVECMULT_THRESHOLD );  // TODO: rethink
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand   sm_;  //!< Sparse matrix of the reduction expression.
   Operation op_;  //!< The reduction operation.
   //**********************************************************************************************

   //**Assignment to vectors***********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a row-wise row-major sparse matrix reduction operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a row-wise row-major
   // sparse matrix reduction expression to a vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto assign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      assign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to vectors**************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a row-wise row-major sparse matrix reduction operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto addAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      addAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to vectors***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a row-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto subAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      subAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a row-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto multAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      multAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to vectors**************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a row-wise row-major sparse matrix reduction operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto divAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      divAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to vectors*******************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a row-wise row-major sparse matrix reduction operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a row-wise row-major
   // sparse matrix reduction expression to a vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to vectors**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a row-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpAddAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpAddAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a row-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpSubAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpSubAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a row-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // row-wise row-major sparse matrix reduction expression to a vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpMultAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpMultAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to vectors**********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a row-wise row-major sparse matrix reduction operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side reduction expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a row-wise
   // row-major sparse matrix reduction expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto smpDivAssign( Vector<VT1,false>& lhs, const SMatReduceExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.sm_ );  // Evaluation of the sparse matrix operand
      smpDivAssign( ~lhs, reduce<1UL>( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Performs a custom reduction operation on the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the given sparse matrix \a sm by means of the
// given reduction operation \a op:

   \code
   blaze::CompressedMatrix<double> a;
   // ... Resizing and initialization
   const double sum = reduce( a, Add() );
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified. Thus the
// behavior is non-deterministic if \a op is not associative or not commutative. Also, the
// operation is undefined if the given reduction operation modifies the values.
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the reduction operation
inline decltype(auto) reduce( const SparseMatrix<MT,SO>& sm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using CT = CompositeType_t<MT>;
   using ET = ElementType_t<MT>;

   const size_t M( (~sm).rows()    );
   const size_t N( (~sm).columns() );

   if( M == 0UL || N == 0UL ) return ET{};

   const size_t iend( SO ? N : M );

   CT tmp( ~sm );

   BLAZE_INTERNAL_ASSERT( tmp.rows()    == M, "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( tmp.columns() == N, "Invalid number of columns" );

   ET redux0{};

   for( size_t i=0UL; i<iend; ++i )
   {
      const auto end( tmp.end(i) );
      auto element( tmp.begin(i) );

      if( element == end ) continue;

      ET redux1( element->value() );
      ++element;

      for( ; element!=end; ++element ) {
         redux1 = op( redux1, element->value() );
      }

      if( i == 0UL ) redux0 = redux1;
      else           redux0 = op( redux0, redux1 );
   }

   return redux0;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation for custom reduction operations on row-major sparse matrices.
// \ingroup sparse_matrix
//
// \param sm The given row-major sparse matrix for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
*/
template< size_t RF      // Reduction flag
        , typename MT    // Type of the sparse matrix
        , typename OP >  // Type of the reduction operation
inline const SMatReduceExpr<RF,MT,OP> reduce_backend( const SparseMatrix<MT,false>& sm, OP op )
{
   return SMatReduceExpr<RF,MT,OP>( ~sm, op );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation for custom reduction operations on column-major sparse matrices.
// \ingroup sparse_matrix
//
// \param sm The given column-major sparse matrix for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
*/
template< size_t RF      // Reduction flag
        , typename MT    // Type of the sparse matrix
        , typename OP >  // Type of the reduction operation
inline decltype(auto) reduce_backend( const SparseMatrix<MT,true>& sm, OP op )
{
   return trans( reduce<1UL-RF>( trans( ~sm ), op ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs a custom reduction operation on the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the rows or columns of the given sparse matrix \a sm by means of the
// given reduction operation \a op. In case the reduction flag \a RF is set to 0, the elements
// of the matrix are reduced along the 0th dimension of the matrix (i.e. column-wise) and the
// result is a row vector. In case \a RF is set to 1, the elements of the matrix are reduced
// along the 1st dimension of the matrix (i.e. row-wise) and the result is a column vector:

   \code
   blaze::CompressedMatrix<double> A;
   // ... Resizing and initialization
   const blaze::DynamicVector<double,rowVector> colsum( reduce<0UL>( A, Add() ) );
   \endcode

   \code
   blaze::CompressedMatrix<double> A;
   // ... Resizing and initialization
   const blaze::DynamicVector<double,columnVector> rowsum( reduce<1UL>( A, Add() ) );
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified. Thus the
// behavior is non-deterministic if \a op is not associative or not commutative. Also, the
// operation is undefined if the given reduction operation modifies the values.
*/
template< size_t RF      // Reduction flag
        , typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the reduction operation
inline decltype(auto) reduce( const SparseMatrix<MT,SO>& sm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_STATIC_ASSERT_MSG( RF < 2UL, "Invalid reduction flag" );

   return reduce_backend<RF>( ~sm, op );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given sparse matrix by means of addition.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the given sparse matrix \a sm by means of
// addition:

   \code
   blaze::CompressedMatrix<int> a{ { 1, 2 }, { 3, 4 } };
   // ... Resizing and initialization
   const int sum = sum( a );  // Results in 10
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sum( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~sm, Add() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given sparse matrix by means of addition.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the rows or columns of the given sparse matrix
// \a sm by means of addition. In case the reduction flag \a RF is set to 0, the elements of the
// matrix are reduced along the 0th dimension of the matrix (i.e. column-wise) and the result is
// a row vector. In case \a RF is set to 1, the elements of the matrix are reduced along the 1st
// dimension of the matrix (i.e. row-wise) and the result is a column vector:

   \code
   blaze::CompressedMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   // ... Resizing and initialization
   const blaze::DynamicVector<int,rowVector> colsum = sum<0UL>( A );  // Results in { 2, 3, 6 }
   \endcode

   \code
   blaze::CompressedMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   // ... Resizing and initialization
   const blaze::DynamicVector<int,columnVector> rowsum = sum<1UL>( A );  // Results in { 3, 8 }
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< size_t RF    // Reduction flag
        , typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sum( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return reduce<RF>( ~sm, Add() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given sparse matrix by means of multiplication.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the given sparse matrix \a sm by means of
// multiplication:

   \code
   blaze::CompressedMatrix<int> a{ { 1, 2 }, { 3, 4 } };
   // ... Resizing and initialization
   const int prod = prod( a );  // Results in 24
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) prod( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~sm, Mult() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given sparse matrix by means of multiplication.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the rows or columns of the given sparse matrix
// \a sm by means of multiplication. In case the reduction flag \a RF is set to 0, the elements
// of the matrix are reduced along the 0th dimension of the matrix (i.e. column-wise) and the
// result is a row vector. In case \a RF is set to 1, the elements of the matrix are reduced
// along the 1st dimension of the matrix (i.e. row-wise) and the result is a column vector:

   \code
   blaze::CompressedMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   // ... Resizing and initialization
   const blaze::DynamicVector<int,rowVector> colsum = sum<0UL>( A );  // Results in { 1, 3, 8 }
   \endcode

   \code
   blaze::CompressedMatrix<int> A{ { 1, 0, 2 }, { 1, 3, 4 } };
   // ... Resizing and initialization
   const blaze::DynamicVector<int,columnVector> rowsum = sum<1UL>( A );  // Results in { 2, 12 }
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< size_t RF    // Reduction flag
        , typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) prod( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return reduce<RF>( ~sm, Mult() );
}
//*************************************************************************************************

} // namespace blaze

#endif
