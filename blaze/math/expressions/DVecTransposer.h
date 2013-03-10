//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecTransposer.h
//  \brief Header file for the dense vector transposer
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECTRANSPOSER_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECTRANSPOSER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECTRANSPOSER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the transposition of a dense vector.
// \ingroup dense_vector_expression
//
// The DVecTransposer class is a wrapper object for the temporary transposition of a dense vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
class DVecTransposer : public DenseVector< DVecTransposer<VT,TF>, TF >
{
 public:
   //**Type definitions****************************************************************************
   typedef DVecTransposer<VT,TF>        This;            //!< Type of this DVecTransposer instance.
   typedef typename VT::TransposeType   ResultType;      //!< Result type for expression template evaluations.
   typedef typename VT::ResultType      TransposeType;   //!< Transpose type for expression template evaluations.
   typedef typename VT::ElementType     ElementType;     //!< Resulting element type.
   typedef typename VT::ReturnType      ReturnType;      //!< Return type for expression template evaluations.
   typedef const This&                  CompositeType;   //!< Data type for composite expression templates.
   typedef typename VT::Reference       Reference;       //!< Reference to a non-constant matrix value.
   typedef typename VT::ConstReference  ConstReference;  //!< Reference to a constant matrix value.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for intrinsic optimization.
   /*! The \a vectorizable compilation flag indicates whether expressions the vector is involved
       in can be optimized via intrinsics. In case the dense vector operand is vectorizable, the
       \a vectorizable compilation flag is set to \a true, otherwise it is set to \a false. */
   enum { vectorizable = VT::vectorizable };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecTransposer class.
   //
   // \param dv The dense vector operand.
   */
   explicit inline DVecTransposer( VT& dv )
      : dv_( dv )  // The dense vector operand
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed value.
   */
   inline Reference operator[]( size_t index ) {
      BLAZE_USER_ASSERT( index < dv_.size(), "Invalid vector access index" );
      return dv_[index];
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the vector elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline ElementType* data() {
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

   //**Reset function******************************************************************************
   /*!\brief Resets the vector elements.
   //
   // \return void
   */
   inline void reset() {
      return dv_.reset();
   }
   //**********************************************************************************************

   //**IsAliased function**************************************************************************
   /*!\brief Returns whether the vector is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the alias corresponds to this vector, \a false if not.
   */
   template< typename Other >  // Data type of the foreign expression
   inline bool isAliased( const Other* alias ) const
   {
      return dv_.isAliased( alias );
   }
   //**********************************************************************************************

   //**Transpose assignment of dense vectors*******************************************************
   /*!\brief Implementation of the transpose assignment of a dense vector.
   //
   // \param rhs The right-hand side dense vector to be assigned.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side dense vector
   inline void assign( const DenseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      const size_t n( size() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t end( n & size_t(-2) );

      for( size_t i=0UL; i<end; i+=2UL ) {
         dv_[i    ] = (~rhs)[i    ];
         dv_[i+1UL] = (~rhs)[i+1UL];
      }
      if( end < n )
         dv_[end] = (~rhs)[end];
   }
   //**********************************************************************************************

   //**Transpose assignment of sparse vectors******************************************************
   /*!\brief Implementation of the transpose assignment of a sparse vector.
   //
   // \param rhs The right-hand side sparse vector to be assigned.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side sparse vector
   inline void assign( const SparseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      typedef typename VT2::ConstIterator  ConstIterator;

      for( ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
         dv_[element->index()] = element->value();
   }
   //**********************************************************************************************

   //**Transpose addition assignment of dense vectors**********************************************
   /*!\brief Implementation of the transpose addition assignment of a dense vector.
   //
   // \param rhs The right-hand side dense vector to be added.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side dense vector
   inline void addAssign( const DenseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      const size_t n( size() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t end( n & size_t(-2) );

      for( size_t i=0UL; i<end; i+=2UL ) {
         dv_[i    ] += (~rhs)[i    ];
         dv_[i+1UL] += (~rhs)[i+1UL];
      }
      if( end < n )
         dv_[end] += (~rhs)[end];
   }
   //**********************************************************************************************

   //**Transpose addition assignment of sparse vectors*********************************************
   /*!\brief Implementation of the transpose addition assignment of a sparse vector.
   //
   // \param rhs The right-hand side sparse vector to be added.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side sparse vector
   inline void addAssign( const SparseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      typedef typename VT2::ConstIterator  ConstIterator;

      for( ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
         dv_[element->index()] += element->value();
   }
   //**********************************************************************************************

   //**Transpose subtraction assignment of dense vectors*******************************************
   /*!\brief Implementation of the transpose subtraction assignment of a dense vector.
   //
   // \param rhs The right-hand side dense vector to be subtracted.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side dense vector
   inline void subAssign( const DenseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      const size_t n( size() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t end( n & size_t(-2) );

      for( size_t i=0UL; i<end; i+=2UL ) {
         dv_[i    ] -= (~rhs)[i    ];
         dv_[i+1UL] -= (~rhs)[i+1UL];
      }
      if( end < n )
         dv_[end] -= (~rhs)[end];
   }
   //**********************************************************************************************

   //**Transpose subtraction assignment of sparse vectors******************************************
   /*!\brief Implementation of the transpose subtraction assignment of a sparse vector.
   //
   // \param rhs The right-hand side sparse vector to be subtracted.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side sparse vector
   inline void subAssign( const SparseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      typedef typename VT2::ConstIterator  ConstIterator;

      for( ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
         dv_[element->index()] -= element->value();
   }
   //**********************************************************************************************

   //**Transpose multiplication assignment of dense vectors****************************************
   /*!\brief Implementation of the transpose multiplication assignment of a dense vector.
   //
   // \param rhs The right-hand side dense vector to be multiplied.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side dense vector
   inline void multAssign( const DenseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      const size_t n( size() );

      BLAZE_INTERNAL_ASSERT( ( n - ( n % 2UL ) ) == ( n & size_t(-2) ), "Invalid end calculation" );
      const size_t end( n & size_t(-2) );

      for( size_t i=0UL; i<end; i+=2UL ) {
         dv_[i    ] *= (~rhs)[i    ];
         dv_[i+1UL] *= (~rhs)[i+1UL];
      }
      if( end < n )
         dv_[end] *= (~rhs)[end];
   }
   //**********************************************************************************************

   //**Transpose multiplication assignment of sparse vectors***************************************
   /*!\brief Implementation of the transpose multiplication assignment of a sparse vector.
   //
   // \param rhs The right-hand side sparse vector to be multiplied.
   // \return void
   //
   // This function must \b NOT be called explicitly! It is used internally for the performance
   // optimized evaluation of expression templates. Calling this function explicitly might result
   // in erroneous results and/or in compilation errors. Instead of using this function use the
   // assignment operator.
   */
   template< typename VT2 >  // Type of the right-hand side dense vector
   inline void multAssign( const SparseVector<VT2,TF>& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );

      BLAZE_INTERNAL_ASSERT( dv_.size() == (~rhs).size(), "Invalid vector sizes" );

      typedef typename VT2::ConstIterator  ConstIterator;

      const VT tmp( dv_ );
      const ConstIterator end( (~rhs).end() );

      dv_.reset();

      for( ConstIterator element=(~rhs).begin(); element!=end; ++element )
         dv_[element->index()] = tmp[element->index()] * element->value();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   VT& dv_;  //!< The dense vector operand.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, !TF );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE( VT );
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the dense vector contained in a DVecTransposer.
// \ingroup dense_vector_expression
//
// \param v The dense vector to be resetted.
// \return void
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline void reset( DVecTransposer<VT,TF>& v )
{
   v.reset();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
