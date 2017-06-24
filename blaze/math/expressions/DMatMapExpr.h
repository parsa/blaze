//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatMapExpr.h
//  \brief Header file for the dense matrix map expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATMAPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATMAPEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/CTransExprTrait.h>
#include <blaze/math/traits/DMatMapExprTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/TDMatMapExprTrait.h>
#include <blaze/math/traits/UnaryMapExprTrait.h>
#include <blaze/math/traits/UnaryMapTrait.h>
#include <blaze/math/typetraits/Columns.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/HasMember.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATMAPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense matrix map() function.
// \ingroup dense_matrix_expression
//
// The DMatMapExpr class represents the compile time expression for the evaluation of a custom
// operation on each element of a dense matrix via the map() function.
*/
template< typename MT  // Type of the dense matrix
        , typename OP  // Type of the custom operation
        , bool SO >    // Storage order
class DMatMapExpr : public DenseMatrix< DMatMapExpr<MT,OP,SO>, SO >
                  , private MatMapExpr
                  , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef ResultType_<MT>    RT;  //!< Result type of the dense matrix expression.
   typedef OppositeType_<MT>  OT;  //!< Opposite type of the dense matrix expression.
   typedef ElementType_<MT>   ET;  //!< Element type of the dense matrix expression.
   typedef ReturnType_<MT>    RN;  //!< Return type of the dense matrix expression.

   //! Definition of the HasSIMDEnabled type trait.
   BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( HasSIMDEnabled, simdEnabled );

   //! Definition of the HasLoad type trait.
   BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( HasLoad, load );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the map expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the map expression. In case the given dense matrix
       expression of type \a MT requires an intermediate evaluation, \a useAssign will be
       set to 1 and the map expression will be evaluated via the \a assign function family.
       Otherwise \a useAssign will be set to 0 and the expression will be evaluated via the
       subscript operator. */
   enum : bool { useAssign = RequiresEvaluation<MT>::value };

   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT2 >
   struct UseAssign {
      enum : bool { value = useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The UseSMPAssign struct is a helper struct for the selection of the parallel evaluation
       strategy. In case either the target matrix or the dense matrix operand is not SMP
       assignable and the matrix operand requires an intermediate evaluation, \a value is set
       to 1 and the expression specific evaluation strategy is selected. Otherwise \a value is
       set to 0 and the default strategy is chosen. */
   template< typename MT2 >
   struct UseSMPAssign {
      enum : bool { value = ( !MT2::smpAssignable || !MT::smpAssignable ) && useAssign };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**SIMD support detection**********************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the detection of the SIMD capabilities of the given custom operation.
   struct UseSIMDEnabledFlag {
      enum : bool { value = OP::template simdEnabled<ET>() };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DMatMapExpr<MT,OP,SO>       This;           //!< Type of this DMatMapExpr instance.
   typedef UnaryMapTrait_<RT,OP>       ResultType;     //!< Result type for expression template evaluations.
   typedef OppositeType_<ResultType>   OppositeType;   //!< Result type with opposite storage order for expression template evaluations.
   typedef TransposeType_<ResultType>  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef ElementType_<ResultType>    ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef decltype( std::declval<OP>()( std::declval<RN>() ) )  ReturnType;

   //! Data type for composite expression templates.
   typedef IfTrue_< useAssign, const ResultType, const DMatMapExpr& >  CompositeType;

   //! Composite data type of the dense matrix expression.
   typedef If_< IsExpression<MT>, const MT, const MT& >  Operand;

   //! Data type of the custom unary operation.
   typedef OP  Operation;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense matrix map expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      typedef std::random_access_iterator_tag  IteratorCategory;  //!< The iterator category.
      typedef ElementType                      ValueType;         //!< Type of the underlying elements.
      typedef ElementType*                     PointerType;       //!< Pointer return type.
      typedef ElementType&                     ReferenceType;     //!< Reference return type.
      typedef ptrdiff_t                        DifferenceType;    //!< Difference between two iterators.

      // STL iterator requirements
      typedef IteratorCategory  iterator_category;  //!< The iterator category.
      typedef ValueType         value_type;         //!< Type of the underlying elements.
      typedef PointerType       pointer;            //!< Pointer return type.
      typedef ReferenceType     reference;          //!< Reference return type.
      typedef DifferenceType    difference_type;    //!< Difference between two iterators.

      //! ConstIterator type of the dense matrix expression.
      typedef ConstIterator_<MT>  IteratorType;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param it Iterator to the initial matrix element.
      // \param op The custom unary operation.
      */
      explicit inline ConstIterator( IteratorType it, OP op )
         : it_( it )  // Iterator to the current matrix element
         , op_( op )  // The custom unary operation
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ConstIterator& operator+=( size_t inc ) {
         it_ += inc;
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
         it_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ConstIterator& operator++() {
         ++it_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator++( int ) {
         return ConstIterator( it_++, op_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ConstIterator& operator--() {
         --it_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator--( int ) {
         return ConstIterator( it_--, op_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReturnType operator*() const {
         return op_( *it_ );
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Access to the SIMD elements of the matrix.
      //
      // \return The resulting SIMD element.
      */
      inline auto load() const noexcept {
         return op_.load( it_.load() );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return it_ == rhs.it_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return it_ != rhs.it_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ConstIterator& rhs ) const {
         return it_ < rhs.it_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ConstIterator& rhs ) const {
         return it_ > rhs.it_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ConstIterator& rhs ) const {
         return it_ <= rhs.it_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ConstIterator& rhs ) const {
         return it_ >= rhs.it_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return it_ - rhs.it_;
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
         return ConstIterator( it.it_ + inc, it.op_ );
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
         return ConstIterator( it.it_ + inc, it.op_ );
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
         return ConstIterator( it.it_ - dec, it.op_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType it_;  //!< Iterator to the current matrix element.
      OP           op_;  //!< The custom unary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum : bool { simdEnabled = MT::simdEnabled &&
                               If_< HasSIMDEnabled<OP>, UseSIMDEnabledFlag, HasLoad<OP> >::value };

   //! Compilation switch for the expression template assignment strategy.
   enum : bool { smpAssignable = MT::smpAssignable };
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   enum : size_t { SIMDSIZE = SIMDTrait<ElementType>::size };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatMapExpr class.
   //
   // \param dm The dense matrix operand of the map expression.
   // \param op The custom unary operation.
   */
   explicit inline DMatMapExpr( const MT& dm, OP op ) noexcept
      : dm_( dm )  // Dense matrix of the map expression
      , op_( op )  // The custom unary operation
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < dm_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.columns(), "Invalid column access index" );
      return op_( dm_(i,j) );
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
      if( i >= dm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= dm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Load function*******************************************************************************
   /*!\brief Access to the SIMD elements of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE auto load( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < dm_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.columns(), "Invalid column access index" );
      BLAZE_INTERNAL_ASSERT( !SO || ( i % SIMDSIZE == 0UL ), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( SO  || ( j % SIMDSIZE == 0UL ), "Invalid column access index" );
      return op_.load( dm_.load(i,j) );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator to the first non-zero element of row \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( dm_.begin(i), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator just past the last non-zero element of row \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( dm_.end(i), op_ );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return dm_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return dm_.columns();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the dense matrix operand.
   //
   // \return The dense matrix operand.
   */
   inline Operand operand() const noexcept {
      return dm_;
   }
   //**********************************************************************************************

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the custom operation.
   //
   // \return A copy of the custom operation.
   */
   inline Operation operation() const {
      return op_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return IsExpression<MT>::value && dm_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return dm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return dm_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return dm_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand   dm_;  //!< Dense matrix of the map expression.
   Operation op_;  //!< The custom unary operation.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix map
   // expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target matrix are identical.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline EnableIf_< And< UseAssign<MT2>
                               , IsSame< UnderlyingNumeric<MT>, UnderlyingNumeric<MT2> > > >
      assign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( ~lhs, rhs.dm_ );
      assign( ~lhs, rhs.op_( ~lhs ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix map
   // expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target vector differ.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline EnableIf_< And< UseAssign<MT2>
                               , Not< IsSame< UnderlyingNumeric<MT>, UnderlyingNumeric<MT2> > > > >
      assign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      assign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix map
   // expression to a sparse matrix. Due to the explicit application of the SFINAE principle,
   // this operator can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order or the target sparse matrix
   friend inline EnableIf_< UseAssign<MT2> >
      assign( SparseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      typedef IfTrue_< SO == SO2, RT, OT >  TmpType;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OT, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT2, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( serial( rhs.dm_ ) );
      assign( ~lhs, rhs.op_( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline EnableIf_< UseAssign<MT2> >
      addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      addAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline EnableIf_< UseAssign<MT2> >
      subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      subAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline EnableIf_< UseAssign<MT2> >
      schurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      schurAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to sparse matrices*************************************************
   // No special implementation for the Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix map expression to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix
   // map expression to a row-major dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected and the underlying
   // numeric data type of the operand and the target matrix are identical.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline EnableIf_< And< UseSMPAssign<MT2>
                               , IsSame< UnderlyingNumeric<MT>, UnderlyingNumeric<MT2> > > >
      smpAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAssign( ~lhs, rhs.dm_ );
      smpAssign( ~lhs, rhs.op_( ~lhs ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix map expression to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix map
   // expression to a row-major dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected and the underlying numeric data type
   // of the operand and the target vector differ.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline EnableIf_< And< UseSMPAssign<MT2>
                               , Not< IsSame< UnderlyingNumeric<MT>, UnderlyingNumeric<MT2> > > > >
      smpAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix map
   // expression to a sparse matrix. Due to the explicit application of the SFINAE principle,
   // this operator can only be selected by the compiler in case the expression specific
   // parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order or the target sparse matrix
   friend inline EnableIf_< UseSMPAssign<MT2> >
      smpAssign( SparseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      typedef IfTrue_< SO == SO2, RT, OT >  TmpType;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OT, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT2, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( rhs.dm_ );
      smpAssign( ~lhs, rhs.op_( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline EnableIf_< UseSMPAssign<MT2> >
      smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpAddAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a
   // dense matrix map expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline EnableIf_< UseSMPAssign<MT2> >
      smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpSubAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // dense matrix map expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline EnableIf_< UseSMPAssign<MT2> >
      smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpSchurAssign( ~lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to sparse matrices*********************************************
   // No special implementation for the SMP Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
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
/*!\brief Evaluates the given custom operation on each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a dm.
//
// The \a map() function evaluates the given custom operation on each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a map() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = map( A, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline const DMatMapExpr<MT,OP,SO> map( const DenseMatrix<MT,SO>& dm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,OP,SO>( ~dm, op );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluates the given custom operation on each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a dm.
//
// The \a forEach() function evaluates the given custom operation on each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a forEach() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = forEach( A, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline const DMatMapExpr<MT,OP,SO> forEach( const DenseMatrix<MT,SO>& dm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,OP,SO>( ~dm, op );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a abs() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a abs() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a abs() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = abs( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Abs,SO> abs( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Abs,SO>( ~dm, Abs() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a floor() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a floor() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a floor() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = floor( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Floor,SO> floor( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Floor,SO>( ~dm, Floor() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a ceil() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a ceil() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a ceil() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = ceil( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Ceil,SO> ceil( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Ceil,SO>( ~dm, Ceil() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a trunc() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a trunc() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a trunc() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = trunc( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Trunc,SO> trunc( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Trunc,SO>( ~dm, Trunc() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a round() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a round() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a round() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = round( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Round,SO> round( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Round,SO>( ~dm, Round() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the complex conjugate of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The conjugate complex of each single element of \a dm.
//
// The \a conj function calculates the complex conjugate of each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a conj function:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Conj,SO> conj( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Conj,SO>( ~dm, Conj() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the conjugate transpose matrix of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The conjugate transpose of \a dm.
//
// The \a ctrans function returns an expression representing the conjugate transpose (also called
// adjoint matrix, Hermitian conjugate matrix or transjugate matrix) of the given input matrix
// \a dm.\n
// The following example demonstrates the use of the \a ctrans function:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = ctrans( A );
   \endcode

// Note that the \a ctrans function has the same effect as manually applying the \a conj and
// \a trans function in any order:

   \code
   B = trans( conj( A ) );  // Computing the conjugate transpose matrix
   B = conj( trans( A ) );  // Computing the conjugate transpose matrix
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const CTransExprTrait_<MT> ctrans( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return trans( conj( ~dm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the real part of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The real part of each single element of \a dm.
//
// The \a real function calculates the real part of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a real function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = real( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Real,SO> real( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Real,SO>( ~dm, Real() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the imaginary part of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The imaginary part of each single element of \a dm.
//
// The \a imag function calculates the imaginary part of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a imag function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = imag( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Imag,SO> imag( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Imag,SO>( ~dm, Imag() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the square root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The square root of each single element of \a dm.
//
// The \a sqrt() function computes the square root of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sqrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sqrt( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Sqrt,SO> sqrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Sqrt,SO>( ~dm, Sqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse square root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$(0..\infty)\f$.
// \return The inverse square root of each single element of \a dm.
//
// The \a invsqrt() function computes the inverse square root of each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invsqrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = invsqrt( A );
   \endcode

// \note All elements are expected to be in the range \f$(0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,InvSqrt,SO> invsqrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,InvSqrt,SO>( ~dm, InvSqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cubic root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The cubic root of each single element of \a dm.
//
// The \a cbrt() function computes the cubic root of each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cbrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = cbrt( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Cbrt,SO> cbrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Cbrt,SO>( ~dm, Cbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cubic root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$(0..\infty)\f$.
// \return The inverse cubic root of each single element of \a dm.
//
// The \a invcbrt() function computes the inverse cubic root of each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invcbrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = invcbrt( A );
   \endcode

// \note All elements are expected to be in the range \f$(0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,InvCbrt,SO> invcbrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,InvCbrt,SO>( ~dm, InvCbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Restricts each single element of the dense matrix \a dm to the range \f$[min..max]\f$.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param min The lower delimiter.
// \param max The upper delimiter.
// \return The matrix with restricted elements.
//
// The \a clamp() function restricts each element of the input matrix \a dm to the range
// \f$[min..max]\f$. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a clamp() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = clamp( A, -1.0, 1.0 );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename DT >  // Type of the delimiters
inline const DMatMapExpr<MT,Clamp<DT>,SO>
   clamp( const DenseMatrix<MT,SO>& dm, const DT& min, const DT& max )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Clamp<DT>,SO>( ~dm, Clamp<DT>( min, max ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the exponential value for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param exp The exponent.
// \return The exponential value of each single element of \a dm.
//
// The \a pow() function computes the exponential value for each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a pow() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = pow( A, 4.2 );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename ET >  // Type of the exponent
inline const DMatMapExpr<MT,Pow<ET>,SO> pow( const DenseMatrix<MT,SO>& dm, ET exp )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ET );

   return DMatMapExpr<MT,Pow<ET>,SO>( ~dm, Pow<ET>( exp ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ e^x \f$ for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// The \a exp() function computes \f$ e^x \f$ for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Exp,SO> exp( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Exp,SO>( ~dm, Exp() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 2^x \f$ for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// The \a exp2() function computes \f$ 2^x \f$ for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp2() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp2( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Exp2,SO> exp2( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Exp2,SO>( ~dm, Exp2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 10^x \f$ for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// The \a exp10() function computes \f$ 10^x \f$ for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp10() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp10( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Exp10,SO> exp10( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Exp10,SO>( ~dm, Exp10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of each single element of \a dm.
//
// The \a log() function computes natural logarithm for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Log,SO> log( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Log,SO>( ~dm, Log() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the binary logarithm for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The binary logarithm of each single element of \a dm.
//
// The \a log2() function computes binary logarithm for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log2() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log2( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Log2,SO> log2( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Log2,SO>( ~dm, Log2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the common logarithm for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The common logarithm of each single element of \a dm.
//
// The \a log10() function computes common logarithm for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log10() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log10( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Log10,SO> log10( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Log10,SO>( ~dm, Log10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The sine of each single element of \a dm.
//
// The \a sin() function computes the sine for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sin() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sin( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Sin,SO> sin( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Sin,SO>( ~dm, Sin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse sine of each single element of \a dm.
//
// The \a asin() function computes the inverse sine for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asin() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = asin( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Asin,SO> asin( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Asin,SO>( ~dm, Asin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The hyperbolic sine of each single element of \a dm.
//
// The \a sinh() function computes the hyperbolic sine for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sinh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sinh( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Sinh,SO> sinh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Sinh,SO>( ~dm, Sinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The inverse hyperbolic sine of each single element of \a dm.
//
// The \a asinh() function computes the inverse hyperbolic sine for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asinh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = asinh( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Asinh,SO> asinh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Asinh,SO>( ~dm, Asinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The cosine of each single element of \a dm.
//
// The \a cos() function computes the cosine for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cos() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = cos( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Cos,SO> cos( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Cos,SO>( ~dm, Cos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse cosine of each single element of \a dm.
//
// The \a acos() function computes the inverse cosine for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acos() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = acos( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Acos,SO> acos( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Acos,SO>( ~dm, Acos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The hyperbolic cosine of each single element of \a dm.
//
// The \a cosh() function computes the hyperbolic cosine for each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cosh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = cosh( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Cosh,SO> cosh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Cosh,SO>( ~dm, Cosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[1..\infty)\f$.
// \return The inverse hyperbolic cosine of each single element of \a dm.
//
// The \a acosh() function computes the inverse hyperbolic cosine for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acosh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = acosh( A );
   \endcode

// \note All elements are expected to be in the range \f$[1..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Acosh,SO> acosh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Acosh,SO>( ~dm, Acosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The tangent of each single element of \a dm.
//
// The \a tan() function computes the tangent for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tan() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = tan( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Tan,SO> tan( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Tan,SO>( ~dm, Tan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The inverse tangent of each single element of \a dm.
//
// The \a atan() function computes the inverse tangent for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atan() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = atan( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Atan,SO> atan( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Atan,SO>( ~dm, Atan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The hyperbolic tangent of each single element of \a dm.
//
// The \a tanh() function computes the hyperbolic tangent for each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tanh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = tanh( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Tanh,SO> tanh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Tanh,SO>( ~dm, Tanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse hyperbolic tangent of each single element of \a dm.
//
// The \a atanh() function computes the inverse hyperbolic tangent for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atanh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = atanh( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Atanh,SO> atanh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Atanh,SO>( ~dm, Atanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the error function for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The error function of each single element of \a dm.
//
// The \a erf() function computes the error function for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erf() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = erf( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Erf,SO> erf( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Erf,SO>( ~dm, Erf() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the complementary error function for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The complementary error function of each single element of \a dm.
//
// The \a erfc() function computes the complementary error function for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erfc() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = erfc( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Erfc,SO> erfc( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatMapExpr<MT,Erfc,SO>( ~dm, Erfc() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Absolute value function for absolute value dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The absolute value dense matrix expression.
// \return The absolute value of each single element of \a dm.
//
// This function implements a performance optimized treatment of the absolute value operation
// on a dense matrix absolute value expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Abs,SO>& abs( const DMatMapExpr<MT,Abs,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a floor() function to a dense matrix \a floor() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a floor() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a floor() operation on
// a dense matrix \a floor() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Floor,SO>& floor( const DMatMapExpr<MT,Floor,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a ceil() function to a dense matrix \a ceil() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a ceil() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a ceil() operation on
// a dense matrix \a ceil() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Ceil,SO>& ceil( const DMatMapExpr<MT,Ceil,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a trunc() function to a dense matrix \a trunc() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a trunc() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a trunc() operation on
// a dense matrix \a trunc() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Trunc,SO>& trunc( const DMatMapExpr<MT,Trunc,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a round() function to a dense matrix \a round() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a round() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a round() operation on
// a dense matrix \a round() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Round,SO>& round( const DMatMapExpr<MT,Round,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for complex conjugate dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The complex conjugate dense matrix expression.
// \return The original dense matrix.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a dense matrix complex conjugate expression. It returns an expression representing the
// original dense matrix:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( conj( A ) );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DMatMapExpr<MT,Conj,SO>::Operand conj( const DMatMapExpr<MT,Conj,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for conjugate transpose dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The conjugate transpose dense matrix expression.
// \return The transpose dense matrix.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a dense matrix conjugate transpose expression. It returns an expression representing the
// transpose of the dense matrix:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( ctrans( A ) );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatTransExpr<MT,!SO> conj( const DMatTransExpr<DMatMapExpr<MT,Conj,SO>,!SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return DMatTransExpr<MT,!SO>( dm.operand().operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Real part function for real part dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The real part dense matrix expression.
// \return The real part of each single element of \a dm.
//
// This function implements a performance optimized treatment of the real part operation on
// a dense matrix real part expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline const DMatMapExpr<MT,Real,SO>& real( const DMatMapExpr<MT,Real,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
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
template< typename MT, typename OP, bool SO >
struct Rows< DMatMapExpr<MT,OP,SO> > : public Rows<MT>
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
template< typename MT, typename OP, bool SO >
struct Columns< DMatMapExpr<MT,OP,SO> > : public Columns<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP, bool SO >
struct IsAligned< DMatMapExpr<MT,OP,SO> >
   : public BoolConstant< IsAligned<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP, bool SO >
struct IsPadded< DMatMapExpr<MT,OP,SO> >
   : public BoolConstant< IsPadded<MT>::value >
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
struct IsSymmetric< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Imag,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,InvSqrt,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,InvCbrt,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, typename ET, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Pow<ET>,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Exp,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Exp2,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Exp10,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Log,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Log2,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Log10,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Sin,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Asin,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Sinh,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Asinh,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Cos,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Acos,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Cosh,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Acosh,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Tan,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Atan,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Tanh,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Atanh,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Erf,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
{};

template< typename MT, bool SO >
struct IsSymmetric< DMatMapExpr<MT,Erfc,SO> >
   : public BoolConstant< IsSymmetric<MT>::value >
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
struct IsHermitian< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Imag,SO> >
   : public BoolConstant< IsBuiltin< ElementType_<MT> >::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,InvSqrt,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,InvCbrt,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, typename ET, bool SO >
struct IsHermitian< DMatMapExpr<MT,Pow<ET>,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Exp,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Exp2,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Exp10,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Log,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Log2,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Log10,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Sin,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Asin,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Sinh,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Asinh,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Cos,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Acos,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Cosh,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Acosh,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Tan,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Atan,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Tanh,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Atanh,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Erf,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
{};

template< typename MT, bool SO >
struct IsHermitian< DMatMapExpr<MT,Erfc,SO> >
   : public BoolConstant< IsHermitian<MT>::value >
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
struct IsLower< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Imag,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Sin,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Asin,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Sinh,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Asinh,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Tan,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Atan,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Tanh,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Atanh,SO> >
   : public BoolConstant< IsLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsLower< DMatMapExpr<MT,Erf,SO> >
   : public BoolConstant< IsLower<MT>::value >
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
struct IsUniLower< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniLower< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
{};

template< typename MT, typename ET, bool SO >
struct IsUniLower< DMatMapExpr<MT,Pow<ET>,SO> >
   : public BoolConstant< IsUniLower<MT>::value >
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
struct IsStrictlyLower< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Imag,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Sin,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Asin,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Sinh,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Asinh,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Tan,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Atan,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Tanh,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Atanh,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyLower< DMatMapExpr<MT,Erf,SO> >
   : public BoolConstant< IsStrictlyLower<MT>::value >
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
struct IsUpper< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Imag,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Sin,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Asin,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Sinh,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Asinh,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Tan,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Atan,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Tanh,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Atanh,SO> >
   : public BoolConstant< IsUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUpper< DMatMapExpr<MT,Erf,SO> >
   : public BoolConstant< IsUpper<MT>::value >
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
struct IsUniUpper< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
{};

template< typename MT, typename ET, bool SO >
struct IsUniUpper< DMatMapExpr<MT,Pow<ET>,SO> >
   : public BoolConstant< IsUniUpper<MT>::value >
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
struct IsStrictlyUpper< DMatMapExpr<MT,Abs,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Floor,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Ceil,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Trunc,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Round,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Conj,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Real,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Imag,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Sqrt,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Cbrt,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Sin,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Asin,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Sinh,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Asinh,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Tan,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Atan,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Tanh,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Atanh,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
{};

template< typename MT, bool SO >
struct IsStrictlyUpper< DMatMapExpr<MT,Erf,SO> >
   : public BoolConstant< IsStrictlyUpper<MT>::value >
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
struct DMatMapExprTrait< DMatMapExpr<MT,Abs,false>, Abs >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatMapExpr<MT,Abs,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Abs,true>, Abs >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatMapExpr<MT,Abs,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatMapExpr<MT,Floor,false>, Floor >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatMapExpr<MT,Floor,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Floor,true>, Floor >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatMapExpr<MT,Floor,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatMapExpr<MT,Ceil,false>, Ceil >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatMapExpr<MT,Ceil,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Ceil,true>, Ceil >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatMapExpr<MT,Ceil,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatMapExpr<MT,Trunc,false>, Trunc >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatMapExpr<MT,Trunc,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Trunc,true>, Trunc >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatMapExpr<MT,Trunc,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatMapExpr<MT,Round,false>, Round >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatMapExpr<MT,Round,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Round,true>, Round >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatMapExpr<MT,Round,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatMapExpr<MT,Conj,false>, Conj >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , Operand_< DMatMapExpr<MT,Conj,false> >
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Conj,true>, Conj >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , Operand_< DMatMapExpr<MT,Conj,true> >
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatTransExpr< DMatMapExpr<MT,Conj,true>, false >, Conj >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatTransExpr<MT,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatTransExpr< DMatMapExpr<MT,Conj,false>, true >, Conj >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatTransExpr<MT,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct DMatMapExprTrait< DMatMapExpr<MT,Real,false>, Real >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsRowMajorMatrix<MT> >
                   , DMatMapExpr<MT,Real,false>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct TDMatMapExprTrait< DMatMapExpr<MT,Real,true>, Real >
{
 public:
   //**********************************************************************************************
   using Type = If_< And< IsDenseMatrix<MT>, IsColumnMajorMatrix<MT> >
                   , DMatMapExpr<MT,Real,true>
                   , INVALID_TYPE >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP, bool SO, bool AF >
struct SubmatrixExprTrait< DMatMapExpr<MT,OP,SO>, AF >
{
 public:
   //**********************************************************************************************
   using Type = UnaryMapExprTrait_< SubmatrixExprTrait_<const MT,AF>, OP >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP, bool SO >
struct RowExprTrait< DMatMapExpr<MT,OP,SO> >
{
 public:
   //**********************************************************************************************
   using Type = UnaryMapExprTrait_< RowExprTrait_<const MT>, OP >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP, bool SO >
struct ColumnExprTrait< DMatMapExpr<MT,OP,SO> >
{
 public:
   //**********************************************************************************************
   using Type = UnaryMapExprTrait_< ColumnExprTrait_<const MT>, OP >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
