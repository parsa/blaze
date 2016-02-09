//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatRealExpr.h
//  \brief Header file for the sparse matrix real part expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATREALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATREALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <iterator>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatRealExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/shims/Real.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/RealExprTrait.h>
#include <blaze/math/traits/RealTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SMatRealExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/TSMatRealExprTrait.h>
#include <blaze/math/typetraits/Columns.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATREALEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the sparse matrix real() function.
// \ingroup sparse_matrix_expression
//
// The SMatRealExpr class represents the compile time expression for the calculation of the
// real part of each element of a sparse matrix via the real() function.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
class SMatRealExpr : public SparseMatrix< SMatRealExpr<MT,SO>, SO >
                   , private MatRealExpr
                   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     RT;  //!< Result type of the sparse matrix expression.
   typedef typename MT::ReturnType     RN;  //!< Return type of the sparse matrix expression.
   typedef typename MT::CompositeType  CT;  //!< Composite type of the sparse matrix expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If the matrix operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   enum { returnExpr = !IsTemporary<RN>::value };

   //! Expression return type for the subscript operator.
   typedef typename RealExprTrait<RN>::Type  ExprReturnType;
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the real part expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the real part expression. In case the given sparse
       matrix expression of type \a MT requires an intermediate evaluation, \a useAssign will
       be set to 1 and the real part expression will be evaluated via the \a assign function
       family. Otherwise \a useAssign will be set to 0 and the expression will be evaluated
       via the subscript operator. */
   enum { useAssign = RequiresEvaluation<MT>::value };

   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct UseAssign {
      enum { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The UseSMPAssign struct is a helper struct for the selection of the parallel evaluation
       strategy. In case either the target matrix or the sparse matrix operand is not SMP
       assignable or the matrix operand requires an intermediate evaluation, \a value is set to
       1 and the expression specific evaluation strategy is selected. Otherwise \a value is set
       to 0 and the default strategy is chosen. */
   template< typename MT2 >
   struct UseSMPAssign {
      enum { value = ( !MT2::smpAssignable || !MT::smpAssignable ) && useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SMatRealExpr<MT,SO>                 This;           //!< Type of this SMatRealExpr instance.
   typedef typename RealTrait<RT>::Type        ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::OppositeType   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const SMatRealExpr& >::Type  CompositeType;

   //! Composite data type of the sparse matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  Operand;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the sparse matrix real part expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the sparse matrix expression.
      typedef ValueIndexPair<ElementType>  Element;

      //! Iterator type of the sparse matrix expression.
      typedef typename RemoveReference<Operand>::Type::ConstIterator  IteratorType;

      typedef std::forward_iterator_tag  IteratorCategory;  //!< The iterator category.
      typedef Element                    ValueType;         //!< Type of the underlying pointers.
      typedef ValueType*                 PointerType;       //!< Pointer return type.
      typedef ValueType&                 ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                  DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying pointers.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      */
      inline ConstIterator( IteratorType it )
         : it_( it )  // Iterator over the elements of the sparse matrix expression
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented expression iterator.
      */
      inline ConstIterator& operator++() {
         ++it_;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline const Element operator*() const {
         return Element( real( it_->value() ), it_->index() );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return Reference to the sparse matrix element at the current iterator position.
      */
      inline const ConstIterator* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse element.
      //
      // \return The current value of the sparse element.
      */
      inline ReturnType value() const {
         return real( it_->value() );
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return it_->index();
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return it_ == rhs.it_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return it_ != rhs.it_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two expression iterators.
      //
      // \param rhs The right-hand side expression iterator.
      // \return The number of elements between the two expression iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return it_ - rhs.it_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType it_;  //!< Iterator over the elements of the sparse matrix expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatRealExpr class.
   //
   // \param sm The sparse matrix operand of the real part expression.
   */
   explicit inline SMatRealExpr( const MT& sm )
      : sm_( sm )  // Sparse matrix of the real part expression
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
      BLAZE_INTERNAL_ASSERT( i < sm_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < sm_.columns(), "Invalid column access index" );
      return real( sm_(i,j) );
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid matrix access index.
   */
   inline ReturnType at( size_t i, size_t j ) const {
      if( i >= sm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= sm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator to the first non-zero element of row \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( sm_.begin(i) );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator just past the last non-zero element of row \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( sm_.end(i) );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const {
      return sm_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const {
      return sm_.columns();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return sm_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return sm_.nonZeros(i);
   }
   //**********************************************************************************************

   //**Find function*******************************************************************************
   /*!\brief Searches for a specific matrix element.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the element in case the index is found, end() iterator otherwise.
   */
   inline ConstIterator find( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( sm_.find( i, j ) );
   }
   //**********************************************************************************************

   //**LowerBound function*************************************************************************
   /*!\brief Returns an iterator to the first index not less then the given index.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the first index not less then the given index, end() iterator otherwise.
   */
   inline ConstIterator lowerBound( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( sm_.lowerBound( i, j ) );
   }
   //**********************************************************************************************

   //**UpperBound function*************************************************************************
   /*!\brief Returns an iterator to the first index greater then the given index.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the first index greater then the given index, end() iterator otherwise.
   */
   inline ConstIterator upperBound( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( sm_.upperBound( i, j ) );
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the sparse matrix operand.
   //
   // \return The sparse matrix operand.
   */
   inline Operand operand() const {
      return sm_;
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
      return sm_.canAlias( alias );
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
      return sm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const {
      return sm_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand sm_;  //!< Sparse matrix of the real part expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix \a real expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side \a real expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix
   // \a real expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation and the numeric element type of the operand is
   // complex.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( DenseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      assign( ~lhs, real( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix \a real expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side \a real expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix \a real
   // expression to a sparse matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the numeric element type of the operand is complex.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      assign( SparseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      (~lhs).reserve( tmp.nonZeros() );
      assign( ~lhs, real( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix \a real expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side \a real expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse
   // matrix \a real expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      addAssign( DenseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      addAssign( ~lhs, real( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix \a real expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side \a real expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix \a real expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline typename EnableIf< UseAssign<MT2> >::Type
      subAssign( DenseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      subAssign( ~lhs, real( tmp ) );
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

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse matrix \a real expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side \a real expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix
   // \a real expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the the expression
   // specific parallel evaluation strategy is selected and the numeric element type of the
   // operand is complex.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseSMPAssign<MT2> >::Type
      smpAssign( DenseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpAssign( ~lhs, real( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   // No special implementation for the SMP assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix \a real expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side \a real expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix \a real expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline typename EnableIf< UseSMPAssign<MT2> >::Type
      smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpAddAssign( ~lhs, real( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix \a real expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side \a real expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix \a real expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline typename EnableIf< UseSMPAssign<MT2> >::Type
      smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const SMatRealExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename RT::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpSubAssign( ~lhs, real( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
   BLAZE_CONSTRAINT_MUST_NOT_BE_BUILTIN_TYPE( typename UnderlyingNumeric<MT>::Type );
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
/*!\brief Returns a matrix containing the real parts of each single element of \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The real part of each single element of \a sm.
//
// The \a real function calculates the real part of each element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a real function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = real( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline const typename RealExprTrait<MT>::Type real( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return typename RealExprTrait<MT>::Type( ~sm );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Real part function for real part sparse matrix expressions.
// \ingroup sparse_matrix
//
// \param sm The real part sparse matrix expression.
// \return The real part of each single element of \a sm.
//
// This function implements a performance optimized treatment of the real part operation on
// a sparse matrix real part expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool TF >    // Transpose flag
inline const SMatRealExpr<MT,TF>& real( const SMatRealExpr<MT,TF>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct Rows< SMatRealExpr<MT,SO> > : public Rows<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct Columns< SMatRealExpr<MT,SO> > : public Columns<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsSymmetric< SMatRealExpr<MT,SO> > : public IsTrue< IsSymmetric<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsHermitian< SMatRealExpr<MT,SO> > : public IsTrue< IsHermitian<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsLower< SMatRealExpr<MT,SO> > : public IsTrue< IsLower<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsUniLower< SMatRealExpr<MT,SO> > : public IsTrue< IsUniLower<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsStrictlyLower< SMatRealExpr<MT,SO> > : public IsTrue< IsStrictlyLower<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsUpper< SMatRealExpr<MT,SO> > : public IsTrue< IsUpper<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsUniUpper< SMatRealExpr<MT,SO> > : public IsTrue< IsUniUpper<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsStrictlyUpper< SMatRealExpr<MT,SO> > : public IsTrue< IsStrictlyUpper<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct SMatRealExprTrait< SMatRealExpr<MT,false> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT>::value && IsRowMajorMatrix<MT>::value
                              , SMatRealExpr<MT,false>
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TSMatRealExprTrait< SMatRealExpr<MT,true> >
{
 public:
   //**********************************************************************************************
   typedef typename SelectType< IsSparseMatrix<MT>::value && IsColumnMajorMatrix<MT>::value
                              , SMatRealExpr<MT,true>
                              , INVALID_TYPE >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool AF >
struct SubmatrixExprTrait< SMatRealExpr<MT,SO>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename RealExprTrait< typename SubmatrixExprTrait<const MT,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct RowExprTrait< SMatRealExpr<MT,SO> >
{
 public:
   //**********************************************************************************************
   typedef typename RealExprTrait< typename RowExprTrait<const MT>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct ColumnExprTrait< SMatRealExpr<MT,SO> >
{
 public:
   //**********************************************************************************************
   typedef typename RealExprTrait< typename ColumnExprTrait<const MT>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
