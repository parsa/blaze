//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatDVecMultExpr.h
//  \brief Header file for the transpose dense matrix/dense vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <boost/cast.hpp>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsBlasCompatible.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/system/BLAS.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/constraints/Double.h>
#include <blaze/util/constraints/Float.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATDVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// The TDMatDVecMultExpr class represents the compile time expression for multiplications
// between column-major dense matrices and dense vectors.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side dense vector
class TDMatDVecMultExpr : public DenseVector< TDMatDVecMultExpr<MT,VT>, false >
                        , private Expression
                        , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::ResultType     MRT;  //!< Result type of the left-hand side dense matrix expression.
   typedef typename VT::ResultType     VRT;  //!< Result type of the right-hand side dense vector expression.
   typedef typename MRT::ElementType   MET;  //!< Element type of the left-hand side dense matrix expression.
   typedef typename VRT::ElementType   VET;  //!< Element type of the right-hand side dense vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the left-hand side dense matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   enum { evaluate = IsComputation<MT>::value && !MT::vectorizable &&
                     IsSame<VET,MET>::value && IsBlasCompatible<VET>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a float and the
       single precision kernel can be used, the nested \a value will be set to 1, otherwise it
       will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseSinglePrecisionKernel {
      enum { value = IsFloat<typename T1::ElementType>::value &&
                     IsFloat<typename T2::ElementType>::value &&
                     IsFloat<typename T3::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a double and the
       double precision kernel can be used, the nested \a value will be set to 1, otherwise it
       will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseDoublePrecisionKernel {
      enum { value = IsDouble<typename T1::ElementType>::value &&
                     IsDouble<typename T2::ElementType>::value &&
                     IsDouble<typename T3::ElementType>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a complex<float>
       and the single precision complex kernel can be used, the nested \a value will be set
       to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseSinglePrecisionComplexKernel {
      typedef complex<float>  Type;
      enum { value = IsSame<typename T1::ElementType,Type>::value &&
                     IsSame<typename T2::ElementType,Type>::value &&
                     IsSame<typename T3::ElementType,Type>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a complex<double>
       and the double precision complex kernel can be used, the nested \a value will be set
       to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseDoublePrecisionComplexKernel {
      typedef complex<double>  Type;
      enum { value = IsSame<typename T1::ElementType,Type>::value &&
                     IsSame<typename T2::ElementType,Type>::value &&
                     IsSame<typename T3::ElementType,Type>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case no optimized BLAS kernel can be used, the nested \a value will be set to 1,
       otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseDefaultKernel {
      enum { value = !BLAZE_BLAS_MODE || ( !UseSinglePrecisionKernel<T1,T2,T3>::value &&
                                           !UseDoublePrecisionKernel<T1,T2,T3>::value &&
                                           !UseSinglePrecisionComplexKernel<T1,T2,T3>::value &&
                                           !UseDoublePrecisionComplexKernel<T1,T2,T3>::value ) };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the matrix type and the two involved vector types are suited for a vectorized
       computation of the matrix/vector multiplication, the nested \value will be set to 1,
       otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseVectorizedDefaultKernel {
      enum { value = T1::vectorizable && T2::vectorizable && T3::vectorizable &&
                     IsSame<typename T1::ElementType,typename T2::ElementType>::value &&
                     IsSame<typename T1::ElementType,typename T3::ElementType>::value &&
                     IntrinsicTrait<typename T1::ElementType>::addition &&
                     IntrinsicTrait<typename T1::ElementType>::multiplication };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef TDMatDVecMultExpr<MT,VT>                    This;           //!< Type of this TDMatDVecMultExpr instance.
   typedef typename MultTrait<MRT,VRT>::Type           ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const ElementType                           ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                            CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   typedef typename SelectType< IsExpression<MT>::value, const MT, const MT& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT>::value, const VT, const VT& >::Type  RightOperand;

   //! Type for the assignment of the left-hand side dense matrix operand.
   typedef typename SelectType< evaluate, const MRT, MCT >::Type  LT;

   //! Type for the assignment of the right-hand side dense vector operand.
   typedef typename SelectType< IsComputation<VT>::value, const VRT, VCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = ( !evaluate && IsComputation<MT>::value && !RequiresEvaluation<MT>::value &&
                       CanAlias<MT>::value ) || ( !IsComputation<VT>::value ) };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDMatDVecMultExpr class.
   //
   // \param mat The left-hand side matrix operand of the multiplication expression.
   // \param vec The right-hand side vector operand of the multiplication expression.
   */
   explicit inline TDMatDVecMultExpr( const MT& mat, const VT& vec )
      : mat_( mat )                                         // Left-hand side dense matrix of the multiplication expression
      , vec_( vec )                                         // Right-hand side dense vector of the multiplication expression
      , end_( ( (mat.columns()-1UL) & size_t(-2) ) + 1UL )  // End of the unrolled calculation loop
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

      if( mat_.columns() != 0UL ) {
         res = mat_(index,0UL) * vec_[0UL];
         for( size_t j=1UL; j<end_; j+=2UL ) {
            res += mat_(index,j) * vec_[j] + mat_(index,j+1) * vec_[j+1UL];
         }
         if( end_ < mat_.columns() ) {
            res += mat_(index,end_) * vec_[end_];
         }
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
   /*!\brief Returns the left-hand side transpose dense matrix operand.
   //
   // \return The left-hand side transpose dense matrix operand.
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
      return ( !evaluate && IsComputation<MT>::value && !RequiresEvaluation<MT>::value &&
               CanAlias<MT>::value && mat_.isAliased( alias ) ) ||
             ( !IsComputation<VT>::value && vec_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand vec_;  //!< Right-hand side dense vector of the multiplication expression.
   const size_t end_;  //!< End of the unrolled calculation loop.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-dense vector multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }
      else if( rhs.mat_.columns() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      if( ( IsComputation<MT>::value && !evaluate ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         TDMatDVecMultExpr::selectDefaultAssignKernel( ~lhs, A, x );
      else
         TDMatDVecMultExpr::selectBlasAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default assignment kernel for the transpose dense matrix-dense
   // vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      BLAZE_INTERNAL_ASSERT( ( M - ( M % 2UL ) ) == ( M & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( M & size_t(-2) );

      for( size_t i=0UL; i<M; ++i ) {
         y[i] = x[0UL] * A(i,0UL);
      }
      for( size_t j=1UL; j<N; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            y[i    ] += x[j] * A(i    ,j);
            y[i+1UL] += x[j] * A(i+1UL,j);
         }
         if( iend < M ) {
            y[iend] += x[j] * A(iend,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.spacing() );
      const size_t N( A.columns() );

      size_t i( 0UL );

      for( ; (i+IT::size*8UL) <= M; i+=IT::size*8UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.get(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.get(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.get(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.get(i+IT::size*7UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
         store( &y[i+IT::size*3UL], xmm4 );
         store( &y[i+IT::size*4UL], xmm5 );
         store( &y[i+IT::size*5UL], xmm6 );
         store( &y[i+IT::size*6UL], xmm7 );
         store( &y[i+IT::size*7UL], xmm8 );
      }
      for( ; (i+IT::size*4UL) <= M; i+=IT::size*4UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
         store( &y[i+IT::size*3UL], xmm4 );
      }
      for( ; (i+IT::size*3UL) <= M; i+=IT::size*3UL ) {
         IntrinsicType xmm1, xmm2, xmm3;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
      }
      for( ; (i+IT::size*2UL) <= M; i+=IT::size*2UL ) {
         IntrinsicType xmm1, xmm2;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i         ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size,j) * x1;
         }
         store( &y[i         ], xmm1 );
         store( &y[i+IT::size], xmm2 );
      }
      if( i < M ) {
         IntrinsicType xmm1;
         for( size_t j=0UL; j<N; ++j ) {
            xmm1 = xmm1 + A.get(i,j) * set( x[j] );
         }
         store( &y[i], xmm1 );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (default)********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDefaultKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (single precision)***********************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-dense vector multiplication for
   //        single precision operands (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for single
   // precision operands based on the BLAS cblas_sgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSinglePrecisionKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_sgemv( CblasColMajor, CblasNoTrans, M, N, 1.0F,
                   A.data(), lda, x.data(), 1, 0.0F, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (double precision)***********************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-dense vector multiplication for
   //        double precision operands (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for double
   // precision operands based on the BLAS cblas_dgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDoublePrecisionKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_dgemv( CblasColMajor, CblasNoTrans, M, N, 1.0,
                   A.data(), lda, x.data(), 1, 0.0, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (single precision complex)***************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-dense vector multiplication for
   //        single precision complex operands (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for single
   // precision complex operands based on the BLAS cblas_cgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSinglePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<float> alpha( 1.0F, 0.0F );
      const complex<float> beta ( 0.0F, 0.0F );

      cblas_cgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (double precision complex)***************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-dense vector multiplication for
   //        double precision complex operands (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for double
   // precision complex operands based on the BLAS cblas_zgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDoublePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<double> alpha( 1.0, 0.0 );
      const complex<double> beta ( 0.0, 0.0 );

      cblas_zgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-dense vector multiplication to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   /*!\brief Addition assignment of a transpose dense matrix-dense vector multiplication to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ) {
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      if( ( IsComputation<MT>::value && !evaluate ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         TDMatDVecMultExpr::selectDefaultAddAssignKernel( ~lhs, A, x );
      else
         TDMatDVecMultExpr::selectBlasAddAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default addition assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      BLAZE_INTERNAL_ASSERT( ( M - ( M % 2UL ) ) == ( M & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( M & size_t(-2) );

      for( size_t j=0UL; j<N; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            y[i    ] += x[j] * A(i    ,j);
            y[i+1UL] += x[j] * A(i+1UL,j);
         }
         if( iend < M ) {
            y[iend] += x[j] * A(iend,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.spacing() );
      const size_t N( A.columns() );

      size_t i( 0UL );

      for( ; (i+IT::size*8UL) <= M; i+=IT::size*8UL ) {
         IntrinsicType xmm1( load( &y[i             ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size    ] ) );
         IntrinsicType xmm3( load( &y[i+IT::size*2UL] ) );
         IntrinsicType xmm4( load( &y[i+IT::size*3UL] ) );
         IntrinsicType xmm5( load( &y[i+IT::size*4UL] ) );
         IntrinsicType xmm6( load( &y[i+IT::size*5UL] ) );
         IntrinsicType xmm7( load( &y[i+IT::size*6UL] ) );
         IntrinsicType xmm8( load( &y[i+IT::size*7UL] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.get(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.get(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.get(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.get(i+IT::size*7UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
         store( &y[i+IT::size*3UL], xmm4 );
         store( &y[i+IT::size*4UL], xmm5 );
         store( &y[i+IT::size*5UL], xmm6 );
         store( &y[i+IT::size*6UL], xmm7 );
         store( &y[i+IT::size*7UL], xmm8 );
      }
      for( ; (i+IT::size*4UL) <= M; i+=IT::size*4UL ) {
         IntrinsicType xmm1( load( &y[i             ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size    ] ) );
         IntrinsicType xmm3( load( &y[i+IT::size*2UL] ) );
         IntrinsicType xmm4( load( &y[i+IT::size*3UL] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
         store( &y[i+IT::size*3UL], xmm4 );
      }
      for( ; (i+IT::size*3UL) <= M; i+=IT::size*3UL ) {
         IntrinsicType xmm1( load( &y[i             ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size    ] ) );
         IntrinsicType xmm3( load( &y[i+IT::size*2UL] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
      }
      for( ; (i+IT::size*2UL) <= M; i+=IT::size*2UL ) {
         IntrinsicType xmm1( load( &y[i         ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i         ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size,j) * x1;
         }
         store( &y[i         ], xmm1 );
         store( &y[i+IT::size], xmm2 );
      }
      if( i < M ) {
         IntrinsicType xmm1( load( &y[i] ) );
         for( size_t j=0UL; j<N; ++j ) {
            xmm1 = xmm1 + A.get(i,j) * set( x[j] );
         }
         store( &y[i], xmm1 );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (default)***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDefaultKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (single precision)**************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose matrix-vector multiplication for single
   //        precision operands (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for single
   // precision operands based on the BLAS cblas_sgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSinglePrecisionKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_sgemv( CblasColMajor, CblasNoTrans, M, N, 1.0F,
                   A.data(), lda, x.data(), 1, 1.0F, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (double precision)**************************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose matrix-vector multiplication for double
   //        precision operands (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for double
   // precision operands based on the BLAS cblas_dgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDoublePrecisionKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_dgemv( CblasColMajor, CblasNoTrans, M, N, 1.0,
                   A.data(), lda, x.data(), 1, 1.0, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (single precision complex)******************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose matrix-vector multiplication for single
   //        precision complex operands (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for single
   // precision complex operands based on the BLAS cblas_cgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSinglePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<float> alpha( 1.0F, 0.0F );
      const complex<float> beta ( 1.0F, 0.0F );

      cblas_cgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (double precision complex)******************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose matrix-vector multiplication for double
   //        precision complex operands (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for double
   // precision complex operands based on the BLAS cblas_zgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDoublePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<double> alpha( 1.0, 0.0 );
      const complex<double> beta ( 1.0, 0.0 );

      cblas_zgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-dense vector multiplication to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ) {
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size()     , "Invalid vector size"       );

      if( ( IsComputation<MT>::value && !evaluate ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         TDMatDVecMultExpr::selectDefaultSubAssignKernel( ~lhs, A, x );
      else
         TDMatDVecMultExpr::selectBlasSubAssignKernel( ~lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      BLAZE_INTERNAL_ASSERT( ( M - ( M % 2UL ) ) == ( M & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( M & size_t(-2) );

      for( size_t j=0UL; j<N; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            y[i    ] -= x[j] * A(i    ,j);
            y[i+1UL] -= x[j] * A(i+1UL,j);
         }
         if( iend < M ) {
            y[iend] -= x[j] * A(iend,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // transpose dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2> >::Type
      selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.spacing() );
      const size_t N( A.columns() );

      size_t i( 0UL );

      for( ; (i+IT::size*8UL) <= M; i+=IT::size*8UL ) {
         IntrinsicType xmm1( load( &y[i             ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size    ] ) );
         IntrinsicType xmm3( load( &y[i+IT::size*2UL] ) );
         IntrinsicType xmm4( load( &y[i+IT::size*3UL] ) );
         IntrinsicType xmm5( load( &y[i+IT::size*4UL] ) );
         IntrinsicType xmm6( load( &y[i+IT::size*5UL] ) );
         IntrinsicType xmm7( load( &y[i+IT::size*6UL] ) );
         IntrinsicType xmm8( load( &y[i+IT::size*7UL] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.get(i             ,j) * x1;
            xmm2 = xmm2 - A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 - A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 - A.get(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 - A.get(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 - A.get(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 - A.get(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 - A.get(i+IT::size*7UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
         store( &y[i+IT::size*3UL], xmm4 );
         store( &y[i+IT::size*4UL], xmm5 );
         store( &y[i+IT::size*5UL], xmm6 );
         store( &y[i+IT::size*6UL], xmm7 );
         store( &y[i+IT::size*7UL], xmm8 );
      }
      for( ; (i+IT::size*4UL) <= M; i+=IT::size*4UL ) {
         IntrinsicType xmm1( load( &y[i             ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size    ] ) );
         IntrinsicType xmm3( load( &y[i+IT::size*2UL] ) );
         IntrinsicType xmm4( load( &y[i+IT::size*3UL] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.get(i             ,j) * x1;
            xmm2 = xmm2 - A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 - A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 - A.get(i+IT::size*3UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
         store( &y[i+IT::size*3UL], xmm4 );
      }
      for( ; (i+IT::size*3UL) <= M; i+=IT::size*3UL ) {
         IntrinsicType xmm1( load( &y[i             ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size    ] ) );
         IntrinsicType xmm3( load( &y[i+IT::size*2UL] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.get(i             ,j) * x1;
            xmm2 = xmm2 - A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 - A.get(i+IT::size*2UL,j) * x1;
         }
         store( &y[i             ], xmm1 );
         store( &y[i+IT::size    ], xmm2 );
         store( &y[i+IT::size*2UL], xmm3 );
      }
      for( ; (i+IT::size*2UL) <= M; i+=IT::size*2UL ) {
         IntrinsicType xmm1( load( &y[i         ] ) );
         IntrinsicType xmm2( load( &y[i+IT::size] ) );
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 - A.get(i         ,j) * x1;
            xmm2 = xmm2 - A.get(i+IT::size,j) * x1;
         }
         store( &y[i         ], xmm1 );
         store( &y[i+IT::size], xmm2 );
      }
      if( i < M ) {
         IntrinsicType xmm1( load( &y[i] ) );
         for( size_t j=0UL; j<N; ++j ) {
            xmm1 = xmm1 - A.get(i,j) * set( x[j] );
         }
         store( &y[i], xmm1 );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (default)********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDefaultKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      selectDefaultSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (single precision)***********************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subtraction assignment of a transpose matrix-vector multiplication for
   //        single precision operands (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for single
   // precision operands based on the BLAS cblas_sgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSinglePrecisionKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_sgemv( CblasColMajor, CblasNoTrans, M, N, -1.0F,
                   A.data(), lda, x.data(), 1, 1.0F, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (double precision)***********************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subtraction assignment of a transpose matrix-vector multiplication for
   //        double precision operands (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for double
   // precision operands based on the BLAS cblas_dgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDoublePrecisionKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_dgemv( CblasColMajor, CblasNoTrans, M, N, -1.0,
                   A.data(), lda, x.data(), 1, 1.0, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (single precision complex)***************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subtraction assignment of a transpose matrix-vector multiplication for
   //        single precision complex operands (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for single
   // precision complex operands based on the BLAS cblas_cgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseSinglePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<float> alpha( -1.0F, 0.0F );
      const complex<float> beta (  1.0F, 0.0F );

      cblas_cgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (double precision complex)***************
#if BLAZE_BLAS_MODE
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subtraction assignment of a transpose matrix-vector multiplication for
   //        double precision complex operands (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication for double
   // precision complex operands based on the BLAS cblas_zgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline typename EnableIf< UseDoublePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<double> alpha( -1.0, 0.0 );
      const complex<double> beta (  1.0, 0.0 );

      cblas_zgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a transpose dense matrix-dense vector multiplication to
   //        a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DVECSCALARMULTEXPR SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expression object for scaled transpose dense matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// This specialization of the DVecScalarMultExpr class represents the compile time expression
// for scaled multiplications between a column-major dense matrix and a non-transpose dense
// vector.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT    // Type of the right-hand side dense vector
        , typename ST >  // Type of the side scalar value
class DVecScalarMultExpr< TDMatDVecMultExpr<MT,VT>, ST, false >
   : public DenseVector< DVecScalarMultExpr< TDMatDVecMultExpr<MT,VT>, ST, false >, false >
   , private Expression
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef TDMatDVecMultExpr<MT,VT>    MVM;  //!< Type of the transpose dense matrix-dense vector multiplication expression.
   typedef typename MVM::ResultType    RES;  //!< Result type of the dense matrix-dense vector multiplication expression.
   typedef typename MT::ResultType     MRT;  //!< Result type of the left-hand side dense matrix expression.
   typedef typename VT::ResultType     VRT;  //!< Result type of the right-hand side dense vector expression.
   typedef typename MRT::ElementType   MET;  //!< Element type of the left-hand side dense matrix expression.
   typedef typename VRT::ElementType   VET;  //!< Element type of the right-hand side dense vector expression.
   typedef typename MT::CompositeType  MCT;  //!< Composite type of the left-hand side dense matrix expression.
   typedef typename VT::CompositeType  VCT;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   enum { evaluate = IsComputation<MT>::value && !MT::vectorizable &&
                     IsSame<VET,MET>::value && IsBlasCompatible<VET>::value };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a float, the scalar
       value is not a complex data type, and the single precision kernel can be used, the nested
       \a value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   struct UseSinglePrecisionKernel {
      enum { value = IsFloat<typename T1::ElementType>::value &&
                     IsFloat<typename T2::ElementType>::value &&
                     IsFloat<typename T3::ElementType>::value &&
                     !IsComplex<T4>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a double, the scalar
       value is not a complex data type and the double precision kernel can be used, the nested
       \a value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   struct UseDoublePrecisionKernel {
      enum { value = IsDouble<typename T1::ElementType>::value &&
                     IsDouble<typename T2::ElementType>::value &&
                     IsDouble<typename T3::ElementType>::value &&
                     !IsComplex<T4>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a complex<float>
       and the single precision complex kernel can be used, the nested \a value will be set to
       1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseSinglePrecisionComplexKernel {
      typedef complex<float>  Type;
      enum { value = IsSame<typename T1::ElementType,Type>::value &&
                     IsSame<typename T2::ElementType,Type>::value &&
                     IsSame<typename T3::ElementType,Type>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the data type of the two involved vectors and the matrix is \a complex<double>
       and the double precision complex kernel can be used, the nested \a value will be set to
       1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   struct UseDoublePrecisionComplexKernel {
      typedef complex<double>  Type;
      enum { value = IsSame<typename T1::ElementType,Type>::value &&
                     IsSame<typename T2::ElementType,Type>::value &&
                     IsSame<typename T3::ElementType,Type>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case no optimized BLAS kernel can be used, the nested \a value will be set to 1,
       otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   struct UseDefaultKernel {
      enum { value = !BLAZE_BLAS_MODE || ( !UseSinglePrecisionKernel<T1,T2,T3,T4>::value &&
                                           !UseDoublePrecisionKernel<T1,T2,T3,T4>::value &&
                                           !UseSinglePrecisionComplexKernel<T1,T2,T3>::value &&
                                           !UseDoublePrecisionComplexKernel<T1,T2,T3>::value ) };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! In case the two involved vector types, the matrix type, and the scalar type are suited
       for a vectorized computation of the scaled vector/matrix multiplication, the nested
       \value will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   struct UseVectorizedDefaultKernel {
      enum { value = T1::vectorizable && T2::vectorizable && T3::vectorizable &&
                     IsSame<typename T1::ElementType,typename T2::ElementType>::value &&
                     IsSame<typename T1::ElementType,typename T3::ElementType>::value &&
                     IsSame<typename T1::ElementType,T4>::value &&
                     IntrinsicTrait<typename T1::ElementType>::addition &&
                     IntrinsicTrait<typename T1::ElementType>::multiplication };
   };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DVecScalarMultExpr<MVM,ST,false>            This;           //!< Type of this DVecScalarMultExpr instance.
   typedef typename MultTrait<RES,ST>::Type            ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType          TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType            ElementType;    //!< Resulting element type.
   typedef typename IntrinsicTrait<ElementType>::Type  IntrinsicType;  //!< Resulting intrinsic element type.
   typedef const ElementType                           ReturnType;     //!< Return type for expression template evaluations.
   typedef const ResultType                            CompositeType;  //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   typedef const TDMatDVecMultExpr<MT,VT>  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsNumeric<ElementType>::value, ElementType, ST >::Type  RightOperand;

   //! Type for the assignment of the left-hand side dense matrix operand.
   typedef typename SelectType< evaluate, const MRT, MCT >::Type  LT;

   //! Type for the assignment of the right-hand side dense vector operand.
   typedef typename SelectType< IsComputation<VT>::value, const VRT, VCT >::Type  RT;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation flag for the detection of aliasing effects.
   enum { canAlias = CanAlias<VT>::value };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecScalarMultExpr class.
   //
   // \param vector The left-hand side dense vector of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   explicit inline DVecScalarMultExpr( const MVM& vector, ST scalar )
      : vector_( vector )  // Left-hand side dense vector of the multiplication expression
      , scalar_( scalar )  // Right-hand side scalar of the multiplication expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < vector_.size(), "Invalid vector access index" );
      return vector_[index] * scalar_;
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return vector_.size();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const {
      return vector_;
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
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vector_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ) {
         reset( ~lhs );
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      if( ( IsComputation<MT>::value && !evaluate ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         DVecScalarMultExpr::selectDefaultAssignKernel( ~lhs, A, x, rhs.scalar_ );
      else
         DVecScalarMultExpr::selectBlasAssignKernel( ~lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the default assignment kernel for the scaled transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      BLAZE_INTERNAL_ASSERT( ( M - ( M % 2UL ) ) == ( M & size_t(-2) ), "Invalid end calculation" );
      const size_t iend( M & size_t(-2) );

      for( size_t i=0UL; i<M; ++i ) {
         y[i] = x[0UL] * A(i,0UL);
      }
      for( size_t j=1UL; j<N; ++j ) {
         for( size_t i=0UL; i<iend; i+=2UL ) {
            y[i    ] += x[j] * A(i    ,j);
            y[i+1UL] += x[j] * A(i+1UL,j);
         }
         if( iend < M ) {
            y[iend] += x[j] * A(iend,j);
         }
      }
      for( size_t i=0UL; i<M; ++i ) {
         y[i] *= scalar;
      }
   }
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors**********************************************
   /*!\brief Vectorized default assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.spacing() );
      const size_t N( A.columns() );

      const IntrinsicType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+IT::size*8UL) <= M; i+=IT::size*8UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.get(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.get(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.get(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.get(i+IT::size*7UL,j) * x1;
         }
         store( &y[i             ], xmm1*factor );
         store( &y[i+IT::size    ], xmm2*factor );
         store( &y[i+IT::size*2UL], xmm3*factor );
         store( &y[i+IT::size*3UL], xmm4*factor );
         store( &y[i+IT::size*4UL], xmm5*factor );
         store( &y[i+IT::size*5UL], xmm6*factor );
         store( &y[i+IT::size*6UL], xmm7*factor );
         store( &y[i+IT::size*7UL], xmm8*factor );
      }
      for( ; (i+IT::size*4UL) <= M; i+=IT::size*4UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
         }
         store( &y[i             ], xmm1*factor );
         store( &y[i+IT::size    ], xmm2*factor );
         store( &y[i+IT::size*2UL], xmm3*factor );
         store( &y[i+IT::size*3UL], xmm4*factor );
      }
      for( ; (i+IT::size*3UL) <= M; i+=IT::size*3UL ) {
         IntrinsicType xmm1, xmm2, xmm3;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
         }
         store( &y[i             ], xmm1*factor );
         store( &y[i+IT::size    ], xmm2*factor );
         store( &y[i+IT::size*2UL], xmm3*factor );
      }
      for( ; (i+IT::size*2UL) <= M; i+=IT::size*2UL ) {
         IntrinsicType xmm1, xmm2;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i         ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size,j) * x1;
         }
         store( &y[i         ], xmm1*factor );
         store( &y[i+IT::size], xmm2*factor );
      }
      if( i < M ) {
         IntrinsicType xmm1;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i,j) * x1;
         }
         store( &y[i], xmm1*factor );
      }
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (default)********************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (single precision)***********************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-dense vector multiplication
   //        for single precision operands (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // single precision operands based on the BLAS cblas_sgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseSinglePrecisionKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_sgemv( CblasColMajor, CblasNoTrans, M, N, scalar,
                   A.data(), lda, x.data(), 1, 0.0F, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (double precision)***********************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-dense vector multiplication
   //        for double precision operands (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // double precision operands based on the BLAS cblas_dgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDoublePrecisionKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_dgemv( CblasColMajor, CblasNoTrans, M, N, scalar,
                   A.data(), lda, x.data(), 1, 0.0, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (single precision complex)***************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-dense vector multiplication
   //        for single precision complex operands (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // single precision complex operands based on the BLAS cblas_cgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseSinglePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST2 );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<float> alpha( scalar );
      const complex<float> beta ( 0.0F, 0.0F );

      cblas_cgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (double precision complex)***************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-dense vector multiplication
   //        for double precision complex operands (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // double precision complex operands based on the BLAS cblas_zgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDoublePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST2 );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<double> alpha( scalar );
      const complex<double> beta ( 0.0, 0.0 );

      cblas_zgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-dense vector multiplication to a sparse
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      assign( ~lhs, tmp );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      if( ( IsComputation<MT>::value && !evaluate ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         DVecScalarMultExpr::selectDefaultAddAssignKernel( ~lhs, A, x, rhs.scalar_ );
      else
         DVecScalarMultExpr::selectBlasAddAssignKernel( ~lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      y.addAssign( A * x * scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors*************************************
   /*!\brief Vectorized default addition assignment of a scaled transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the scaled
   // transpose dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.spacing() );
      const size_t N( A.columns() );

      const IntrinsicType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+IT::size*8UL) <= M; i+=IT::size*8UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.get(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.get(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.get(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.get(i+IT::size*7UL,j) * x1;
         }
         store( &y[i             ], load( &y[i             ] ) + xmm1*factor );
         store( &y[i+IT::size    ], load( &y[i+IT::size    ] ) + xmm2*factor );
         store( &y[i+IT::size*2UL], load( &y[i+IT::size*2UL] ) + xmm3*factor );
         store( &y[i+IT::size*3UL], load( &y[i+IT::size*3UL] ) + xmm4*factor );
         store( &y[i+IT::size*4UL], load( &y[i+IT::size*4UL] ) + xmm5*factor );
         store( &y[i+IT::size*5UL], load( &y[i+IT::size*5UL] ) + xmm6*factor );
         store( &y[i+IT::size*6UL], load( &y[i+IT::size*6UL] ) + xmm7*factor );
         store( &y[i+IT::size*7UL], load( &y[i+IT::size*7UL] ) + xmm8*factor );
      }
      for( ; (i+IT::size*4UL) <= M; i+=IT::size*4UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
         }
         store( &y[i             ], load( &y[i             ] ) + xmm1*factor );
         store( &y[i+IT::size    ], load( &y[i+IT::size    ] ) + xmm2*factor );
         store( &y[i+IT::size*2UL], load( &y[i+IT::size*2UL] ) + xmm3*factor );
         store( &y[i+IT::size*3UL], load( &y[i+IT::size*3UL] ) + xmm4*factor );
      }
      for( ; (i+IT::size*3UL) <= M; i+=IT::size*3UL ) {
         IntrinsicType xmm1, xmm2, xmm3;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
         }
         store( &y[i             ], load( &y[i             ] ) + xmm1*factor );
         store( &y[i+IT::size    ], load( &y[i+IT::size    ] ) + xmm2*factor );
         store( &y[i+IT::size*2UL], load( &y[i+IT::size*2UL] ) + xmm3*factor );
      }
      for( ; (i+IT::size*2UL) <= M; i+=IT::size*2UL ) {
         IntrinsicType xmm1, xmm2;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i         ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size,j) * x1;
         }
         store( &y[i         ], load( &y[i         ] ) + xmm1*factor );
         store( &y[i+IT::size], load( &y[i+IT::size] ) + xmm2*factor );
      }
      if( i < M ) {
         IntrinsicType xmm1;
         for( size_t j=0UL; j<N; ++j ) {
            xmm1 = xmm1 + A.get(i,j) * set( x[j] );
         }
         store( &y[i], load( &y[i] ) + xmm1*factor );
      }
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (default)***********************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (single precision)**************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for single precision operands (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // single precision operands based on the BLAS cblas_sgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseSinglePrecisionKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_sgemv( CblasColMajor, CblasNoTrans, M, N, scalar,
                   A.data(), lda, x.data(), 1, 1.0F, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (double precision)**************************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for double precision operands (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // double precision operands based on the BLAS cblas_dgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDoublePrecisionKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_dgemv( CblasColMajor, CblasNoTrans, M, N, scalar,
                   A.data(), lda, x.data(), 1, 1.0, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (single precision complex)******************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for single precision complex operands (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // single precision complex operands based on the BLAS cblas_cgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseSinglePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST2 );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<float> alpha( scalar );
      const complex<float> beta ( 1.0F, 0.0F );

      cblas_cgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (double precision complex)******************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for double precision complex operands (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // double precision complex operands based on the BLAS cblas_zgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDoublePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST2 );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<double> alpha( scalar );
      const complex<double> beta ( 1.0, 0.0 );

      cblas_zgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a scaled
   // dense transpose matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      typename MVM::LeftOperand  left ( rhs.vector_.leftOperand()  );
      typename MVM::RightOperand right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).size() , "Invalid vector size"       );

      if( ( IsComputation<MT>::value && !evaluate ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         DVecScalarMultExpr::selectDefaultSubAssignKernel( ~lhs, A, x, rhs.scalar_ );
      else
         DVecScalarMultExpr::selectBlasSubAssignKernel( ~lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename DisableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      y.subAssign( A * x * scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors**********************************
   /*!\brief Vectorized default subtraction assignment of a scaled transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the scaled
   // transpose dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseVectorizedDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      typedef IntrinsicTrait<ElementType>  IT;

      const size_t M( A.spacing() );
      const size_t N( A.columns() );

      const IntrinsicType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+IT::size*8UL) <= M; i+=IT::size*8UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
            xmm5 = xmm5 + A.get(i+IT::size*4UL,j) * x1;
            xmm6 = xmm6 + A.get(i+IT::size*5UL,j) * x1;
            xmm7 = xmm7 + A.get(i+IT::size*6UL,j) * x1;
            xmm8 = xmm8 + A.get(i+IT::size*7UL,j) * x1;
         }
         store( &y[i             ], load( &y[i             ] ) - xmm1*factor );
         store( &y[i+IT::size    ], load( &y[i+IT::size    ] ) - xmm2*factor );
         store( &y[i+IT::size*2UL], load( &y[i+IT::size*2UL] ) - xmm3*factor );
         store( &y[i+IT::size*3UL], load( &y[i+IT::size*3UL] ) - xmm4*factor );
         store( &y[i+IT::size*4UL], load( &y[i+IT::size*4UL] ) - xmm5*factor );
         store( &y[i+IT::size*5UL], load( &y[i+IT::size*5UL] ) - xmm6*factor );
         store( &y[i+IT::size*6UL], load( &y[i+IT::size*6UL] ) - xmm7*factor );
         store( &y[i+IT::size*7UL], load( &y[i+IT::size*7UL] ) - xmm8*factor );
      }
      for( ; (i+IT::size*4UL) <= M; i+=IT::size*4UL ) {
         IntrinsicType xmm1, xmm2, xmm3, xmm4;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
            xmm4 = xmm4 + A.get(i+IT::size*3UL,j) * x1;
         }
         store( &y[i             ], load( &y[i             ] ) - xmm1*factor );
         store( &y[i+IT::size    ], load( &y[i+IT::size    ] ) - xmm2*factor );
         store( &y[i+IT::size*2UL], load( &y[i+IT::size*2UL] ) - xmm3*factor );
         store( &y[i+IT::size*3UL], load( &y[i+IT::size*3UL] ) - xmm4*factor );
      }
      for( ; (i+IT::size*3UL) <= M; i+=IT::size*3UL ) {
         IntrinsicType xmm1, xmm2, xmm3;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i             ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size    ,j) * x1;
            xmm3 = xmm3 + A.get(i+IT::size*2UL,j) * x1;
         }
         store( &y[i             ], load( &y[i             ] ) - xmm1*factor );
         store( &y[i+IT::size    ], load( &y[i+IT::size    ] ) - xmm2*factor );
         store( &y[i+IT::size*2UL], load( &y[i+IT::size*2UL] ) - xmm3*factor );
      }
      for( ; (i+IT::size*2UL) <= M; i+=IT::size*2UL ) {
         IntrinsicType xmm1, xmm2;
         for( size_t j=0UL; j<N; ++j ) {
            const IntrinsicType x1( set( x[j] ) );
            xmm1 = xmm1 + A.get(i         ,j) * x1;
            xmm2 = xmm2 + A.get(i+IT::size,j) * x1;
         }
         store( &y[i         ], load( &y[i         ] ) - xmm1*factor );
         store( &y[i+IT::size], load( &y[i+IT::size] ) - xmm2*factor );
      }
      if( i < M ) {
         IntrinsicType xmm1;
         for( size_t j=0UL; j<N; ++j ) {
            xmm1 = xmm1 + A.get(i,j) * set( x[j] );
         }
         store( &y[i], load( &y[i] ) - xmm1*factor );
      }
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (default)********************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDefaultKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      selectDefaultSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (single precision)***********************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for single precision operands (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // single precision operands based on the BLAS cblas_sgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseSinglePrecisionKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_sgemv( CblasColMajor, CblasNoTrans, M, N, -scalar,
                   A.data(), lda, x.data(), 1, 1.0F, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (double precision)***********************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for double precision operands (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // double precision operands based on the BLAS cblas_dgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDoublePrecisionKernel<VT1,MT1,VT2,ST2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );

      cblas_dgemv( CblasColMajor, CblasNoTrans, M, N, -scalar,
                   A.data(), lda, x.data(), 1, 1.0, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (single precision complex)***************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for single precision complex operands (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // single precision complex operands based on the BLAS cblas_cgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseSinglePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST2 );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<float> alpha( -scalar );
      const complex<float> beta ( 1.0F, 0.0F );

      cblas_cgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (double precision complex)***************
#if BLAZE_BLAS_MODE
   /*!\brief BLAS-based subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication for double precision complex operands (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication for
   // double precision complex operands based on the BLAS cblas_zgemv() function.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline typename EnableIf< UseDoublePrecisionComplexKernel<VT1,MT1,VT2> >::Type
      selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      using boost::numeric_cast;

      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
      BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST2 );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
      BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

      const int M  ( numeric_cast<int>( A.rows() )    );
      const int N  ( numeric_cast<int>( A.columns() ) );
      const int lda( numeric_cast<int>( A.spacing() ) );
      const complex<double> alpha( -scalar );
      const complex<double> beta ( 1.0, 0.0 );

      cblas_zgemv( CblasColMajor, CblasNoTrans, M, N, &alpha,
                   A.data(), lda, x.data(), 1, &beta, y.data(), 1 );
   }
#endif
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}*=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( ResultType );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( MVM );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( MVM );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ST );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a column-major dense matrix and a dense
//        vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side column-major dense matrix for the multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a column-major dense matrix and a dense
// vector:

   \code
   using blaze::columnMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<double,columnMajor> A;
   blaze::DynamicVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved element types \a T1::ElementType and \a T2::ElementType. Both the
// dense matrix type \a T1 and the dense vector type \a T2 as well as the two element types
// \a T1::ElementType and \a T2::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side dense vector
inline const typename DisableIf< IsMatMatMultExpr<T1>, TDMatDVecMultExpr<T1,T2> >::Type
   operator*( const DenseMatrix<T1,true>& mat, const DenseVector<T2,false>& vec )
{
   if( (~mat).columns() != (~vec).size() )
      throw std::invalid_argument( "Matrix and vector sizes do not match" );

   return TDMatDVecMultExpr<T1,T2>( ~mat, ~vec );
}
//*************************************************************************************************

} // namespace blaze

#endif
