//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecTransExpr.h
//  \brief Header file for the dense vector transpose expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECTRANSEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECTRANSEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/DVecTransposer.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SVecTransposer.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EmptyType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECTRANSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector transpositions.
// \ingroup dense_vector_expression
//
// The DVecTransExpr class represents the compile time expression for transpositions of
// dense vectors.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DVecTransExpr : public DenseVector< DVecTransExpr<VT,TF>, TF >
                    , private Expression
                    , private SelectType< IsComputation<VT>::value, Computation, EmptyType >::Type
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::CompositeType  CT;  //!< Composite type of the dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the transposition expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the transposition expression. In case the given dense vector
       expression of type \a VT requires an intermediate evaluation, \a useAssign will be set
       to 1 and the transposition expression will be evaluated via the \a assign function
       family. Otherwise \a useAssign will be set to 0 and the expression will be evaluated
       via the subscript operator. */
   enum { useAssign = !IsReference<CT>::value };
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
   typedef DVecTransExpr<VT,TF>        This;           //!< Type of this DVecTransExpr instance.
   typedef typename VT::TransposeType  ResultType;     //!< Result type for expression template evaluations.
   typedef typename VT::ResultType     TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename VT::ElementType    ElementType;    //!< Resulting element type.
   typedef typename VT::ReturnType     ReturnType;     //!< Return type for expression template evaluations.

   //! Resulting intrinsic element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;

   //! Data type for composite expression templates.
   typedef typename SelectType< useAssign, const ResultType, const DVecTransExpr& >::Type  CompositeType;

   //! Composite data type of the dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  Operand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = VT::vectorizable };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<VT>::value };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecTransExpr class.
   //
   // \param dv The dense vector operand of the transposition expression.
   */
   explicit inline DVecTransExpr( const VT& dv )
      : dv_( dv )  // Dense vector of the transposition expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < dv_.size(), "Invalid vector access index" );
      return dv_[index];
   }
   //**********************************************************************************************

   //**Get function********************************************************************************
   /*!\brief Access to the intrinsic elements of the vector.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   inline IntrinsicType get( size_t index ) const {
      typedef IntrinsicTrait<ElementType>  IT;
      BLAZE_INTERNAL_ASSERT( index < dv_.size()      , "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % IT::size == 0UL , "Invalid vector access index" );
      return dv_.get( index );
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the vector elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline const ElementType* data() const {
      return dv_.data();
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return dv_.size();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the dense vector operand.
   //
   // \return The dense vector operand.
   */
   inline Operand operand() const {
      return dv_;
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
      return dv_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand dv_;  //!< Dense vector of the transposition expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector transposition expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side transposition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector
   // transposition expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( DenseVector<VT2,TF>& lhs, const DVecTransExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      DVecTransposer<VT2,!TF> tmp( ~lhs );
      assign( tmp, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector transposition expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side transposition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector
   // transposition expression to a sparse vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      assign( SparseVector<VT2,TF>& lhs, const DVecTransExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      SVecTransposer<VT2,!TF> tmp( ~lhs );
      assign( tmp, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector transposition expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side transposition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // vector transposition expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this operator can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      addAssign( DenseVector<VT2,TF>& lhs, const DVecTransExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      DVecTransposer<VT2,!TF> tmp( ~lhs );
      addAssign( tmp, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector transposition expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side transposition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector transposition expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      subAssign( DenseVector<VT2,TF>& lhs, const DVecTransExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      DVecTransposer<VT2,!TF> tmp( ~lhs );
      subAssign( tmp, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector transposition expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side transposition expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector transposition expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline typename EnableIf< UseAssign<VT2> >::Type
      multAssign( DenseVector<VT2,TF>& lhs, const DVecTransExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      DVecTransposer<VT2,!TF> tmp( ~lhs );
      multAssign( tmp, rhs.dv_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Trans function******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Calculating the transpose of a transpose dense vector.
   // \ingroup dense_vector
   //
   // \param dv The dense vector to be (re-)transposed.
   // \return The transpose of the transpose vector.
   //
   // This function returns an expression representing the transpose of a transpose dense vector:

      \code
      using blaze::columnVector;

      blaze::DynamicVector<double,columnVector> a, b;
      // ... Resizing and initialization
      b = trans( trans( a ) );
      \endcode
   */
   template< typename VT2  // Type of the dense vector
           , bool TF2 >    // Transpose flag of the dense vector
   friend inline Operand trans( const DVecTransExpr<VT2,TF2>& dv )
   {
      return dv.dv_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, !TF );
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
/*!\brief Calculation of the transpose of the given dense vector.
// \ingroup dense_vector
//
// \param dv The dense vector to be transposed.
// \return The transpose of the dense vector.
//
// This function returns an expression representing the transpose of the given dense vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<double,columnVector> a;
   blaze::DynamicVector<double,rowVector>    b;
   // ... Resizing and initialization
   b = trans( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline const DVecTransExpr<VT,!TF> trans( const DenseVector<VT,TF>& dv )
{
   return DVecTransExpr<VT,!TF>( ~dv );
}
//*************************************************************************************************

} // namespace blaze

#endif
