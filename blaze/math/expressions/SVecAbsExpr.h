//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecAbsExpr.h
//  \brief Header file for the sparse vector absolute value expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECABSEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECABSEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <iterator>
#include <boost/type_traits/remove_reference.hpp>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AbsExprTrait.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SVECABSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the sparse vector abs() function.
// \ingroup sparse_vector_expression
//
// The SVecAbsExpr class represents the compile time expression for the calculation of the
// absolute value of each element of a sparse vector via the abs() function.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
class SVecAbsExpr : public SparseVector< SVecAbsExpr<VT,TF>, TF >
                  , private Expression
                  , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::ResultType     RT;  //!< Result type of the sparse vector expression.
   typedef typename VT::ReturnType     RN;  //!< Return type of the sparse vector expression.
   typedef typename VT::CompositeType  CT;  //!< Composite type of the sparse vector expression.
   typedef typename VT::TransposeType  TT;  //!< Transpose type of the sparse vector expression.
   typedef typename VT::ElementType    ET;  //!< Element type of the sparse vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If the vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   enum { returnExpr = !IsTemporary<RN>::value };

   //! Expression return type for the subscript operator.
   typedef typename AbsExprTrait<RN>::Type  ExprReturnType;
   //**********************************************************************************************

   //**Evaluation strategy*************************************************************************
   //! Compilation switch for the evaluation strategy of the absolute value expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the absolute value expression. In case the given sparse vector
       expression of type \a VT requires an intermediate evaluation, \a useAssign will be set
       to 1 and the absolute value expression will be evaluated via the \a assign function
       family. Otherwise \a useAssign will be set to 0 and the expression will be evaluated
       via the subscript operator. */
   enum { useAssign = !IsReference<CT>::value };

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
   typedef SVecAbsExpr<VT,TF>  This;           //!< Type of this SVecAbsExpr instance.
   typedef RT                  ResultType;     //!< Result type for expression template evaluations.
   typedef TT                  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef ET                  ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const SVecAbsExpr& >::Type  CompositeType;

   //! Composite data type of the sparse vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  Operand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<VT>::value };
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the sparse vector absolute value expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the sparse vector expression.
      typedef SparseElement<ElementType>  Element;

      //! Iterator type of the sparse vector expression.
      typedef typename boost::remove_reference<Operand>::type::ConstIterator  IteratorType;

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
         : it_( it )  // Iterator over the elements of the sparse vector expression
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
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline const Element operator*() const {
         using std::abs;
         return Element( abs( it_->value() ), it_->index() );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
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
         using std::abs;
         return abs( it_->value() );
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
      IteratorType it_;  //!< Iterator over the elements of the sparse vector expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SVecAbsExpr class.
   //
   // \param sv The sparse vector operand of the absolute value expression.
   */
   explicit inline SVecAbsExpr( const VT& sv )
      : sv_( sv )  // Sparse vector of the absolute value expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      using std::abs;
      BLAZE_INTERNAL_ASSERT( index < sv_.size(), "Invalid vector access index" );
      return abs( sv_[index] );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the sparse vector.
   //
   // \return Iterator to the first non-zero element of the sparse vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( sv_.begin() );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the sparse vector.
   //
   // \return Iterator just past the last non-zero element of the sparse vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( sv_.end() );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return sv_.size();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse vector.
   //
   // \return The number of non-zero elements in the sparse vector.
   */
   inline size_t nonZeros() const {
      return sv_.nonZeros();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the sparse vector operand.
   //
   // \return The sparse vector operand.
   */
   inline Operand operand() const {
      return sv_;
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
      return sv_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand sv_;  //!< Sparse vector of the absolute value expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector abs expression to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side abs expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector
   // abs expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( DenseVector<VT2,TF>& lhs, const SVecAbsExpr& rhs )
   {
      using std::abs;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( ~lhs, rhs.sv_ );

      const size_t size( rhs.size() );
      for( size_t i=0UL; i<size; ++i ) {
         (~lhs)[i] = abs( (~lhs)[i] );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector abs expression to a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side abs expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector
   // abs expression to a sparse vector. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( SparseVector<VT2,TF>& lhs, const SVecAbsExpr& rhs )
   {
      using std::abs;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typedef typename VT2::Iterator  Iterator;

      assign( ~lhs, rhs.sv_ );

      const Iterator end( (~lhs).end() );
      for( Iterator element=(~lhs).begin(); element!=end; ++element ) {
         element->value() = abs( element->value() );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse vector abs expression to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side abs expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse
   // vector abs expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      addAssign( DenseVector<VT2,TF>& lhs, const SVecAbsExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      addAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse vector abs expression to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side abs expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // vector abs expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      subAssign( DenseVector<VT2,TF>& lhs, const SVecAbsExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      subAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse vector abs expression to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side abs expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // vector abs expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      multAssign( DenseVector<VT2,TF>& lhs, const SVecAbsExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a vector containing the absolute values of each single element of \a sv.
// \ingroup sparse_vector
//
// \param sv The integral sparse input vector.
// \return The absolute value of each single element of \a sv.
//
// The \a abs function calculates the absolute value of each element of the sparse input vector
// \a sv. The operator returns an expression representing this operation.\n
// The following example demonstrates the use of the \a abs function:

   \code
   blaze::CompressedVector<double> a, b;
   // ... Resizing and initialization
   b = abs( a );
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline const SVecAbsExpr<VT,TF> abs( const SparseVector<VT,TF>& sv )
{
   return SVecAbsExpr<VT,TF>( ~sv );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Absolute value operator for absolute value sparse vector expressions.
// \ingroup sparse_vector
//
// \param sv The absolute value sparse vector expression.
// \return The absolute value of each single element of \a sv.
//
// This operator implements a performance optimized treatment of the absolute value operation
// on a sparse vector absolute value expression.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline const SVecAbsExpr<VT,TF>& abs( const SVecAbsExpr<VT,TF>& sv )
{
   return sv;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
