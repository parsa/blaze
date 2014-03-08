//=================================================================================================
/*!
//  \file blaze/math/dense/StaticMatrix.h
//  \brief Header file for the implementation of a fixed-size matrix
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

#ifndef _BLAZE_MATH_DENSE_STATICMATRIX_H_
#define _BLAZE_MATH_DENSE_STATICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <blaze/math/dense/DenseIterator.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MathTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/util/AlignedArray.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/IsVectorizable.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup static_matrix StaticMatrix
// \ingroup dense_matrix
*/
/*!\brief Efficient implementation of a fixed-sized matrix.
// \ingroup static_matrix
//
// The StaticMatrix class template is the representation of a fixed-size matrix with statically
// allocated elements of arbitrary type. The type of the elements, the number of rows and columns
// and the storage order of the matrix can be specified via the four template parameters:

   \code
   template< typename Type, size_t M, size_t N, bool SO >
   class StaticMatrix;
   \endcode

//  - Type: specifies the type of the matrix elements. StaticMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - M   : specifies the total number of rows of the matrix.
//  - N   : specifies the total number of columns of the matrix. Note that it is expected
//          that StaticMatrix is only used for tiny and small matrices.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::rowMajor.
//
// Depending on the storage order, the matrix elements are either stored in a row-wise fashion
// or in a column-wise fashion. Given the 2x3 matrix

                          \f[\left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          4 & 5 & 6 \\
                          \end{array}\right)\f]\n

// in case of row-major order the elements are stored in the order

                          \f[\left(\begin{array}{*{6}{c}}
                          1 & 2 & 3 & 4 & 5 & 6. \\
                          \end{array}\right)\f]

// In case of column-major order the elements are stored in the order

                          \f[\left(\begin{array}{*{6}{c}}
                          1 & 4 & 2 & 5 & 3 & 6. \\
                          \end{array}\right)\f]

// The use of StaticMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combination of row-major and
// column-major dense and sparse matrices with fitting element types. The following example gives
// an impression of the use of StaticMatrix:

   \code
   using blaze::StaticMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   StaticMatrix<double,2UL,3UL,rowMajor> A;   // Default constructed, non-initialized, row-major 2x3 matrix
   A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;  // Initialization of the first row
   A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;  // Initialization of the second row

   DynamicMatrix<float,2UL,3UL,columnMajor> B;  // Default constructed column-major single precision 2x3 matrix
   B(0,0) = 1.0; B(0,1) = 3.0; B(0,2) = 5.0;    // Initialization of the first row
   B(1,0) = 2.0; B(1,1) = 4.0; B(1,2) = 6.0;    // Initialization of the second row

   CompressedMatrix<float>     C( 2, 3 );  // Empty row-major sparse single precision matrix
   StaticMatrix<float,3UL,2UL> D( 4.0F );  // Directly, homogeneously initialized single precision 3x2 matrix

   StaticMatrix<double,2UL,3UL,rowMajor>    E( A );  // Creation of a new row-major matrix as a copy of A
   StaticMatrix<double,2UL,2UL,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;     // Matrix addition and assignment to a row-major matrix
   E = A - C;     // Matrix subtraction and assignment to a column-major matrix
   F = A * D;     // Matrix multiplication between two matrices of different element types

   A *= 2.0;      // In-place scaling of matrix A
   E  = 2.0 * B;  // Scaling of matrix B
   E  = B * 2.0;  // Scaling of matrix B

   E += A - B;    // Addition assignment
   E -= A + C;    // Subtraction assignment
   F *= A * D;    // Multiplication assignment
   \endcode

// In order to provide maximum performance StaticMatrix is guaranteed to be properly aligned in
// memory based on the alignment restrictions of the specified element type. Note however that
// this enforced alignment can cause problems when StaticMatrix is used in conjunction with a
// std::vector or a similar container:

   \code
   typedef blaze::StaticMatrix<double,3UL,3UL>  MatrixType;

   std::vector< MatrixType > vec( 10 );  // Potentially causes segmentation faults
   \endcode

// The core problem is that the allocated dynamic memory of the container potentially does not
// satisfy the alignment restrictions of the StaticMatrix. For that reason, the AlignedAllocator
// can be used:

   \code
   typedef blaze::StaticMatrix<double,3UL,3UL>  MatrixType;
   typedef blaze::AlignedAllocator<MatrixType>  AllocatorType;

   std::vector< MatrixType, AllocatorType > vec( 10 );  // Guaranteed to work
   \endcode
*/
template< typename Type                    // Data type of the matrix
        , size_t M                         // Number of rows
        , size_t N                         // Number of columns
        , bool SO = defaultStorageOrder >  // Storage order
class StaticMatrix : public DenseMatrix< StaticMatrix<Type,M,N,SO>, SO >
{
 private:
   //**Type definitions****************************************************************************
   typedef IntrinsicTrait<Type>  IT;  //!< Intrinsic trait for the matrix element type.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Alignment adjustment
   enum { NN = N + ( IT::size - ( N % IT::size ) ) % IT::size };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef StaticMatrix<Type,M,N,SO>   This;            //!< Type of this StaticMatrix instance.
   typedef This                        ResultType;      //!< Result type for expression template evaluations.
   typedef StaticMatrix<Type,M,N,!SO>  OppositeType;    //!< Result type with opposite storage order for expression template evaluations.
   typedef StaticMatrix<Type,N,M,!SO>  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef Type                        ElementType;     //!< Type of the matrix elements.
   typedef typename IT::Type           IntrinsicType;   //!< Intrinsic type of the matrix elements.
   typedef const Type&                 ReturnType;      //!< Return type for expression template evaluations.
   typedef const This&                 CompositeType;   //!< Data type for composite expression templates.
   typedef Type&                       Reference;       //!< Reference to a non-constant matrix value.
   typedef const Type&                 ConstReference;  //!< Reference to a constant matrix value.
   typedef DenseIterator<Type>         Iterator;        //!< Iterator over non-constant elements.
   typedef DenseIterator<const Type>   ConstIterator;   //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for intrinsic optimization.
   /*! The \a vectorizable compilation flag indicates whether expressions the matrix is involved
       in can be optimized via intrinsics. In case the element type of the matrix is a vectorizable
       data type, the \a vectorizable compilation flag is set to \a true, otherwise it is set to
       \a false. */
   enum { vectorizable = IsVectorizable<Type>::value };

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline StaticMatrix();
   explicit inline StaticMatrix( const Type& init );

   template< typename Other > explicit inline StaticMatrix( size_t m, size_t n, const Other* array );
   template< typename Other > explicit inline StaticMatrix( const Other (&array)[M][N] );

                                        inline StaticMatrix( const StaticMatrix& m );
   template< typename Other, bool SO2 > inline StaticMatrix( const StaticMatrix<Other,M,N,SO2>& m );
   template< typename MT   , bool SO2 > inline StaticMatrix( const Matrix<MT,SO2>& m );

   inline StaticMatrix( const Type& v1, const Type& v2 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7, const Type& v8 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7, const Type& v8, const Type& v9 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7, const Type& v8, const Type& v9, const Type& v10 );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Type*          data  ();
   inline const Type*    data  () const;
   inline Type*          data  ( size_t i );
   inline const Type*    data  ( size_t i ) const;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   template< typename Other >
   inline StaticMatrix& operator=( const Other (&array)[M][N] );

                                        inline StaticMatrix& operator= ( const Type& set );
                                        inline StaticMatrix& operator= ( const StaticMatrix& rhs );
   template< typename Other, bool SO2 > inline StaticMatrix& operator= ( const StaticMatrix<Other,M,N,SO2>& rhs );
   template< typename MT   , bool SO2 > inline StaticMatrix& operator= ( const Matrix<MT,SO2>& rhs );
   template< typename MT   , bool SO2 > inline StaticMatrix& operator+=( const Matrix<MT,SO2>& rhs );
   template< typename MT   , bool SO2 > inline StaticMatrix& operator-=( const Matrix<MT,SO2>& rhs );
   template< typename MT   , bool SO2 > inline StaticMatrix& operator*=( const Matrix<MT,SO2>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, StaticMatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, StaticMatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t        rows() const;
                              inline size_t        columns() const;
                              inline size_t        spacing() const;
                              inline size_t        capacity() const;
                              inline size_t        capacity( size_t i ) const;
                              inline size_t        nonZeros() const;
                              inline size_t        nonZeros( size_t i ) const;
                              inline void          reset();
                              inline void          reset( size_t i );
                              inline StaticMatrix& transpose();
   template< typename Other > inline StaticMatrix& scale( const Other& scalar );
                              inline void          swap( StaticMatrix& m ) /* throw() */;
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IsRowMajorMatrix<MT>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IntrinsicTrait<Type>::addition &&
                     IsRowMajorMatrix<MT>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IntrinsicTrait<Type>::subtraction &&
                     IsRowMajorMatrix<MT>::value };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool isAligned() const;

   inline IntrinsicType load  ( size_t i, size_t j ) const;
   inline IntrinsicType loadu ( size_t i, size_t j ) const;
   inline void          store ( size_t i, size_t j, const IntrinsicType& value );
   inline void          storeu( size_t i, size_t j, const IntrinsicType& value );
   inline void          stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT, bool SO2 >
   inline typename DisableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,SO2>& rhs );

   template< typename MT, bool SO2 >
   inline typename EnableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,SO2>& rhs );

   template< typename MT > inline void assign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT, bool SO2 >
   inline typename DisableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,SO2>& rhs );

   template< typename MT, bool SO2 >
   inline typename EnableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,SO2>& rhs );

   template< typename MT > inline void addAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT, bool SO2 >
   inline typename DisableIf< VectorizedSubAssign<MT> >::Type
      subAssign( const DenseMatrix<MT,SO2>& rhs );

   template< typename MT, bool SO2 >
   inline typename EnableIf< VectorizedSubAssign<MT> >::Type
      subAssign( const DenseMatrix<MT,SO2>& rhs );

   template< typename MT > inline void subAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,!SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   AlignedArray<Type,M*NN> v_;  //!< The statically allocated matrix elements.
                                /*!< Access to the matrix elements is gained via the function call
                                     operator. In case of row-major order the memory layout of the
                                     elements is
                                     \f[\left(\begin{array}{*{5}{c}}
                                     0            & 1             & 2             & \cdots & N-1         \\
                                     N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                                     \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                                     M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
                                     \end{array}\right)\f]. */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_STATIC_ASSERT( NN % IT::size == 0UL );
   BLAZE_STATIC_ASSERT( NN >= N );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for StaticMatrix.
//
// All matrix elements are initialized to the default value (i.e. 0 for integral data types).
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix()
{
   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M*NN; ++i )
         v_[i] = Type();
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogenous initialization of all elements.
//
// \param init Initial value for all matrix elements.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& init )
{
   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j )
         v_[i*NN+j] = init;

      if( IsNumeric<Type>::value ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param array Dynamic array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a dynamic array:

   \code
   using blaze::rowMajor;

   int* array = new int[6];
   // ... Initialization of the dynamic array
   blaze::StaticMatrix<int,3,4,rowMajor> v( array, 2UL, 3UL );
   delete[] array;
   \endcode

// The matrix is initialized with the values from the given array. Missing values are initialized
// with default values. In case the specified number of rows and/or columns exceeds the maximum
// number of rows/column of the static matrix (i.e. \a m is larger than M or \a n is larger than
// N), a \a std::invalid_argument exception is thrown.\n
// Note that it is expected that the given \a array has at least \a m by \a n elements. Providing
// an array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the initialization array
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( size_t m, size_t n, const Other* array )
{
   if( m > M || n > N )
      throw std::invalid_argument( "Invalid setup of static matrix" );

   for( size_t i=0UL; i<m; ++i ) {
      for( size_t j=0UL; j<n; ++j )
         v_[i*NN+j] = array[i*n+j];

      if( IsNumeric<Type>::value ) {
         for( size_t j=n; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }

   for( size_t i=m; i<M; ++i ) {
      for( size_t j=0UL; j<NN; ++j )
         v_[i*NN+j] = Type();
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all matrix elements.
//
// \param array \f$ M \times N \f$ dimensional array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a static array:

   \code
   using blaze::rowMajor;

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::StaticMatrix<int,3,3,rowMajor> A( init );
   \endcode

// The matrix is initialized with the values from the given array. Missing values are initialized
// with default values (as e.g. the value 6 in the example).
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the initialization array
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Other (&array)[M][N] )
{
   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j )
         v_[i*NN+j] = array[i][j];

      if( IsNumeric<Type>::value ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for StaticMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined in order to enable/facilitate NRV optimization.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const StaticMatrix& m )
{
   for( size_t i=0UL; i<M*NN; ++i )
      v_[i] = m.v_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different StaticMatrix instances.
//
// \param m Matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , size_t M        // Number of rows
        , size_t N        // Number of columns
        , bool SO >       // Storage order
template< typename Other  // Data type of the foreign matrix
        , bool SO2 >      // Storage order of the foreign matrix
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const StaticMatrix<Other,M,N,SO2>& m )
{
   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j )
         v_[i*NN+j] = m(i,j);

      if( IsNumeric<Type>::value ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
// \exception std::invalid_argument Invalid setup of static matrix.
//
// This constructor initializes the static matrix from the given matrix. In case the size of
// the given matrix does not match the size of the static matrix (i.e. the number of rows is
// not M or the number of columns is not N), a \a std::invalid_argument exception is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the foreign matrix
        , bool SO2 >     // Storage order of the foreign matrix
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Matrix<MT,SO2>& m )
{
   using blaze::assign;

   if( (~m).rows() != M || (~m).columns() != N )
      throw std::invalid_argument( "Invalid setup of static matrix" );

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=( IsSparseMatrix<MT>::value )?( 0UL ):( N ); j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }

   assign( *this, ~m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 2 \f$ and \f$ 2 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 2 \f$
// and \f$ 2 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{2}{c}}
                     1 & 2 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,2,false> A( 1, 2 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2 )
{
   BLAZE_STATIC_ASSERT( M*N == 2UL );

   // Initialization of a 1x2 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
   }

   // Initialization of a 2x1 matrix
   else {
      v_[0UL] = v1;
      v_[ NN] = v2;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 3 \f$ and \f$ 3 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 3 \f$
// and \f$ 3 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{3}{c}}
                     1 & 2 & 3 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,3,false> A( 1, 2, 3 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3 )
{
   BLAZE_STATIC_ASSERT( M*N == 3UL );

   // Initialization of a 1x3 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
   }

   // Initialization of a 3x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 4 \f$, \f$ 2 \times 2 \f$, and \f$ 4 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 4 \f$,
// \f$ 2 \times 2 \f$ and \f$ 3 \times 1 \f$ matrix. The following example demonstrates this by
// creating the matrix

                     \f[\left(\begin{array}{*{2}{c}}
                     1 & 2 \\
                     3 & 4 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,2,false> A( 1, 2, 3, 4 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2,
                                                const Type& v3, const Type& v4 )
{
   BLAZE_STATIC_ASSERT( M*N == 4UL );

   // Initialization of a 1x4 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
   }

   // Initialization of a 2x2 matrix
   else if( M == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[NN    ] = v3;
      v_[NN+1UL] = v4;
   }

   // Initialization of a 4x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 5 \f$ and \f$ 5 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 5 \f$,
// and \f$ 5 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{5}{c}}
                     1 & 2 & 3 & 4 & 5 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,5,false> A( 1, 2, 3, 4, 5 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                const Type& v4, const Type& v5 )
{
   BLAZE_STATIC_ASSERT( M*N == 5UL );

   // Initialization of a 1x5 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
   }

   // Initialization of a 5x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
      v_[4UL*NN] = v5;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 6 \f$, \f$ 2 \times 3 \f$, \f$ 3 \times 2 \f$, and
//        \f$ 6 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 6 \f$,
// \f$ 2 \times 3 \f$, \f$ 3 \times 2 \f$, and \f$ 6 \times 1 \f$ matrix. The following example
// demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{3}{c}}
                     1 & 2 & 3 \\
                     4 & 5 & 6 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,3,false> A( 1, 2, 3, 4, 5, 6 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                const Type& v4, const Type& v5, const Type& v6 )
{
   BLAZE_STATIC_ASSERT( M*N == 6UL );

   // Initialization of a 1x6 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
   }

   // Initialization of a 2x3 matrix
   else if( M == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[   2UL] = v3;
      v_[NN    ] = v4;
      v_[NN+1UL] = v5;
      v_[NN+2UL] = v6;
   }

   // Initialization of a 3x2 matrix
   else if( M == 3UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[    NN    ] = v3;
      v_[    NN+1UL] = v4;
      v_[2UL*NN    ] = v5;
      v_[2UL*NN+1UL] = v6;
   }

   // Initialization of a 6x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
      v_[4UL*NN] = v5;
      v_[5UL*NN] = v6;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 7 \f$ and \f$ 7 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 7 \f$,
// and \f$ 7 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{7}{c}}
                     1 & 2 & 3 & 4 & 5 & 6 & 7 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,7,false> A( 1, 2, 3, 4, 5, 6, 7 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                const Type& v4, const Type& v5, const Type& v6,
                                                const Type& v7 )
{
   BLAZE_STATIC_ASSERT( M*N == 7UL );

   // Initialization of a 1x7 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
   }

   // Initialization of a 7x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
      v_[4UL*NN] = v5;
      v_[5UL*NN] = v6;
      v_[6UL*NN] = v7;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 8 \f$ and \f$ 8 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
// \param v8 The initializer for the eigth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 8 \f$,
// \f$ 2 \times 4 \f$, \f$ 4 \times 2 \f$, and \f$ 8 \times 1 \f$ matrix. The following example
// demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{4}{c}}
                     1 & 2 & 3 & 4 \\
                     5 & 6 & 7 & 8 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,4,false> A( 1, 2, 3, 4, 5, 6, 7, 8 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                const Type& v4, const Type& v5, const Type& v6,
                                                const Type& v7, const Type& v8 )
{
   BLAZE_STATIC_ASSERT( M*N == 8UL );

   // Initialization of a 1x8 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
      v_[7UL] = v8;
   }

   // Initialization of a 2x4 matrix
   else if( M == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[   2UL] = v3;
      v_[   3UL] = v4;
      v_[NN    ] = v5;
      v_[NN+1UL] = v6;
      v_[NN+2UL] = v7;
      v_[NN+3UL] = v8;
   }

   // Initialization of a 4x2 matrix
   else if( M == 4UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[    NN    ] = v3;
      v_[    NN+1UL] = v4;
      v_[2UL*NN    ] = v5;
      v_[2UL*NN+1UL] = v6;
      v_[3UL*NN    ] = v7;
      v_[3UL*NN+1UL] = v8;
   }

   // Initialization of a 8x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
      v_[4UL*NN] = v5;
      v_[5UL*NN] = v6;
      v_[6UL*NN] = v7;
      v_[7UL*NN] = v8;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 9 \f$, \f$ 3 \times 3 \f$, and \f$ 9 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
// \param v8 The initializer for the eigth matrix element.
// \param v9 The initializer for the ninth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 9 \f$,
// \f$ 3 \times 3 \f$, and \f$ 9 \times 1 \f$ matrix. The following example demonstrates this by
// creating the matrix

                     \f[\left(\begin{array}{*{3}{c}}
                     1 & 2 & 3 \\
                     4 & 5 & 6 \\
                     7 & 8 & 9 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,3,3,false> A( 1, 2, 3, 4, 5, 6, 7, 8, 9 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                const Type& v4, const Type& v5, const Type& v6,
                                                const Type& v7, const Type& v8, const Type& v9 )
{
   BLAZE_STATIC_ASSERT( M*N == 9UL );

   // Initialization of a 1x9 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
      v_[7UL] = v8;
      v_[8UL] = v9;
   }

   // Initialization of a 3x3 matrix
   else if( M == 3UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[       2UL] = v3;
      v_[    NN    ] = v4;
      v_[    NN+1UL] = v5;
      v_[    NN+2UL] = v6;
      v_[2UL*NN    ] = v7;
      v_[2UL*NN+1UL] = v8;
      v_[2UL*NN+2UL] = v9;
   }

   // Initialization of a 9x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
      v_[4UL*NN] = v5;
      v_[5UL*NN] = v6;
      v_[6UL*NN] = v7;
      v_[7UL*NN] = v8;
      v_[8UL*NN] = v9;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for \f$ 1 \times 10 \f$, \f$ 2 \times 5 \f$, \f$ 5 \times 2 \f$, and
//        \f$ 10 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
// \param v8 The initializer for the eigth matrix element.
// \param v9 The initializer for the ninth matrix element.
// \param v10 The initializer for the tenth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 10 \f$,
// \f$ 2 \times 5 \f$, \f$ 5 \times 2 \f$, and \f$ 10 \times 1 \f$ matrix. The following example
// demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{5}{c}}
                     1 & 2 & 3 & 4 & 5  \\
                     6 & 7 & 8 & 8 & 10 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,5,false> A( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                const Type& v4, const Type& v5, const Type& v6,
                                                const Type& v7, const Type& v8, const Type& v9,
                                                const Type& v10 )
{
   BLAZE_STATIC_ASSERT( M*N == 10UL );

   // Initialization of a 1x10 matrix
   if( M == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
      v_[7UL] = v8;
      v_[8UL] = v9;
      v_[9UL] = v10;
   }

   // Initialization of a 2x5 matrix
   else if( M == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[   2UL] = v3;
      v_[   3UL] = v4;
      v_[   4UL] = v5;
      v_[NN    ] = v6;
      v_[NN+1UL] = v7;
      v_[NN+2UL] = v8;
      v_[NN+3UL] = v9;
      v_[NN+4UL] = v10;
   }

   // Initialization of a 5x2 matrix
   else if( M == 5UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[    NN    ] = v3;
      v_[    NN+1UL] = v4;
      v_[2UL*NN    ] = v5;
      v_[2UL*NN+1UL] = v6;
      v_[3UL*NN    ] = v7;
      v_[3UL*NN+1UL] = v8;
      v_[4UL*NN    ] = v9;
      v_[4UL*NN+1UL] = v10;
   }

   // Initialization of a 10x1 matrix
   else {
      v_[   0UL] = v1;
      v_[    NN] = v2;
      v_[2UL*NN] = v3;
      v_[3UL*NN] = v4;
      v_[4UL*NN] = v5;
      v_[5UL*NN] = v6;
      v_[6UL*NN] = v7;
      v_[7UL*NN] = v8;
      v_[8UL*NN] = v9;
      v_[9UL*NN] = v10;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=N; j<NN; ++j )
            v_[i*NN+j] = Type();
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::Reference
   StaticMatrix<Type,M,N,SO>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i<M, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<N, "Invalid column access index" );
   return v_[i*NN+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return Reference-to-const to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::ConstReference
   StaticMatrix<Type,M,N,SO>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i<M, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<N, "Invalid column access index" );
   return v_[i*NN+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the static matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The static matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline Type* StaticMatrix<Type,M,N,SO>::data()
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the static matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The static matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline const Type* StaticMatrix<Type,M,N,SO>::data() const
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements of row/column \a i.
//
// \param i The row/column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row/column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline Type* StaticMatrix<Type,M,N,SO>::data( size_t i )
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return v_ + i*NN;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements of row/column \a i.
//
// \param i The row/column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row/column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline const Type* StaticMatrix<Type,M,N,SO>::data( size_t i ) const
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return v_ + i*NN;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// This function returns a row/column iterator to the first element of row/column \a i. In case
// the storage order is set to \a rowMajor the function returns an iterator to the first element
// of row \a i, in case the storage flag is set to \a columnMajor the function returns an iterator
// to the first element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::Iterator
   StaticMatrix<Type,M,N,SO>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return Iterator( v_ + i*NN );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// This function returns a row/column iterator to the first element of row/column \a i. In case
// the storage order is set to \a rowMajor the function returns an iterator to the first element
// of row \a i, in case the storage flag is set to \a columnMajor the function returns an iterator
// to the first element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::ConstIterator
   StaticMatrix<Type,M,N,SO>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*NN );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// This function returns a row/column iterator to the first element of row/column \a i. In case
// the storage order is set to \a rowMajor the function returns an iterator to the first element
// of row \a i, in case the storage flag is set to \a columnMajor the function returns an iterator
// to the first element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::ConstIterator
   StaticMatrix<Type,M,N,SO>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*NN );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// This function returns an row/column iterator just past the last element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator just past
// the last element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator just past the last element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::Iterator
   StaticMatrix<Type,M,N,SO>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return Iterator( v_ + i*NN + N );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// This function returns an row/column iterator just past the last element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator just past
// the last element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator just past the last element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::ConstIterator
   StaticMatrix<Type,M,N,SO>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*NN + N );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// This function returns an row/column iterator just past the last element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator just past
// the last element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator just past the last element of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::ConstIterator
   StaticMatrix<Type,M,N,SO>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < M, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*NN + N );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Array assignment to all matrix elements.
//
// \param array \f$ M \times N \f$ dimensional array for the assignment.
// \return Reference to the assigned matrix.
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::rowMajor;

   const real init[3][3] = { { 1, 2, 3 },
                             { 4, 5 },
                             { 7, 8, 9 } };
   blaze::StaticMatrix<int,3UL,3UL,rowMajor> A;
   A = init;
   \endcode

// The matrix is assigned the values from the given array. Missing values are initialized with
// default values (as e.g. the value 6 in the example).
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the initialization array
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator=( const Other (&array)[M][N] )
{
   for( size_t i=0UL; i<M; ++i )
      for( size_t j=0UL; j<N; ++j )
         v_[i*NN+j] = array[i][j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Homogenous assignment to all matrix elements.
//
// \param set Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator=( const Type& set )
{
   for( size_t i=0UL; i<M; ++i )
      for( size_t j=0UL; j<N; ++j )
         v_[i*NN+j] = set;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for StaticMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// Explicit definition of a copy assignment operator for performance reasons.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator=( const StaticMatrix& rhs )
{
   using blaze::assign;

   assign( *this, ~rhs );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different StaticMatrix instances.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
*/
template< typename Type   // Data type of the matrix
        , size_t M        // Number of rows
        , size_t N        // Number of columns
        , bool SO >       // Storage order
template< typename Other  // Data type of the foreign matrix
        , bool SO2 >      // Storage order of the foreign matrix
inline StaticMatrix<Type,M,N,SO>&
   StaticMatrix<Type,M,N,SO>::operator=( const StaticMatrix<Other,M,N,SO2>& rhs )
{
   using blaze::assign;

   assign( *this, ~rhs );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid assignment to static matrix.
//
// This constructor initializes the matrix as a copy of the given matrix. In case the
// number of rows of the given matrix is not M or the number of columns is not N, a
// \a std::invalid_argument exception is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator=( const Matrix<MT,SO2>& rhs )
{
   using blaze::assign;

   if( (~rhs).rows() != M || (~rhs).columns() != N )
      throw std::invalid_argument( "Invalid assignment to static matrix" );

   if( (~rhs).canAlias( this ) ) {
      StaticMatrix tmp( ~rhs );
      swap( tmp );
   }
   else {
      if( IsSparseMatrix<MT>::value )
         reset();
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator+=( const Matrix<MT,SO2>& rhs )
{
   using blaze::addAssign;

   if( (~rhs).rows() != M || (~rhs).columns() != N )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      StaticMatrix tmp( ~rhs );
      addAssign( *this, tmp );
   }
   else {
      addAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator-=( const Matrix<MT,SO2>& rhs )
{
   using blaze::subAssign;

   if( (~rhs).rows() != M || (~rhs).columns() != N )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      StaticMatrix tmp( ~rhs );
      subAssign( *this, tmp );
   }
   else {
      subAssign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::operator*=( const Matrix<MT,SO2>& rhs )
{
   if( M != N || (~rhs).rows() != M || (~rhs).columns() != M )
      throw std::invalid_argument( "Matrix sizes do not match" );

   StaticMatrix tmp( *this * (~rhs) );
   return this->operator=( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a matrix and
//        a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, StaticMatrix<Type,M,N,SO> >::Type&
   StaticMatrix<Type,M,N,SO>::operator*=( Other rhs )
{
   using blaze::assign;

   assign( *this, (*this) * rhs );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a matrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the matrix.
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, StaticMatrix<Type,M,N,SO> >::Type&
   StaticMatrix<Type,M,N,SO>::operator/=( Other rhs )
{
   using blaze::assign;

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   assign( *this, (*this) / rhs );
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the matrix.
//
// \return The number of rows of the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::rows() const
{
   return M;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::columns() const
{
   return N;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows.
//
// \return The spacing between the beginning of two rows.
//
// This function returns the spacing between the beginning of two rows, i.e. the total number
// of elements of a row.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::spacing() const
{
   return NN;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::capacity() const
{
   return M*NN;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row/column.
//
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// This function returns the current capacity of the specified row/column. In case the
// storage order is set to \a rowMajor the function returns the capacity of row \a i,
// in case the storage flag is set to \a columnMajor the function returns the capacity
// of column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::capacity( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return NN;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the total number of non-zero elements in the matrix
//
// \return The number of non-zero elements in the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<M; ++i )
      for( size_t j=0UL; j<N; ++j )
         if( !isDefault( v_[i*NN+j] ) )
            ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column.
// In case the storage order is set to \a rowMajor the function returns the number of non-zero
// elements in row \a i, in case the storage flag is set to \a columnMajor the function returns
// the number of non-zero elements in column \a i.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline size_t StaticMatrix<Type,M,N,SO>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( (i+1UL)*NN );
   size_t nonzeros( 0UL );

   for( size_t j=i*NN; j<jend; ++j )
      if( !isDefault( v_[j] ) )
         ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void StaticMatrix<Type,M,N,SO>::reset()
{
   using blaze::reset;

   for( size_t i=0UL; i<M; ++i )
      for( size_t j=0UL; j<N; ++j )
         reset( v_[i*NN+j] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column to the default initial values.
//
// \param i The index of the row/column.
// \return void
//
// This function resets the values in the specified row/column to their default value. In case
// the storage order is set to \a rowMajor the function resets the values in row \a i, in case
// the storage order is set to \a columnMajor the function resets the values in column \a i.
// Note that the capacity of the row/column remains unchanged.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void StaticMatrix<Type,M,N,SO>::reset( size_t i )
{
   using blaze::reset;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   for( size_t j=0UL; j<N; ++j )
      reset( v_[i*NN+j] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transposing the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::transpose()
{
   using std::swap;

   for( size_t i=1UL; i<M; ++i )
      for( size_t j=0UL; j<i; ++j )
         swap( v_[i*NN+j], v_[j*NN+i] );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A*=s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the matrix.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline StaticMatrix<Type,M,N,SO>& StaticMatrix<Type,M,N,SO>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<M; ++i )
      for( size_t j=0UL; j<N; ++j )
         v_[i*NN+j] *= scalar;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two static matrices.
//
// \param m The matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void StaticMatrix<Type,M,N,SO>::swap( StaticMatrix& m ) /* throw() */
{
   using std::swap;

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         swap( v_[i*NN+j], m(i,j) );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the matrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address can alias with the matrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool StaticMatrix<Type,M,N,SO>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address is aliased with the matrix. In contrast
// to the conAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool StaticMatrix<Type,M,N,SO>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix is properly aligned in memory.
//
// \return \a true in case the matrix is aligned, \a false if not.
//
// This function returns whether the matrix is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of each row/column of the matrix are guaranteed to conform
// to the alignment restrictions of the element type \a Type.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline bool StaticMatrix<Type,M,N,SO>::isAligned() const
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::IntrinsicType
   StaticMatrix<Type,M,N,SO>::load( size_t i, size_t j ) const
{
   using blaze::load;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= NN , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   return load( &v_[i*NN+j] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline typename StaticMatrix<Type,M,N,SO>::IntrinsicType
   StaticMatrix<Type,M,N,SO>::loadu( size_t i, size_t j ) const
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= NN , "Invalid column access index" );

   return loadu( &v_[i*NN+j] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void StaticMatrix<Type,M,N,SO>::store( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::store;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= NN , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   store( &v_[i*NN+j], value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void StaticMatrix<Type,M,N,SO>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= NN, "Invalid column access index" );

   storeu( &v_[i*NN+j], value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific intrinsic element of the
// dense matrix. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the column index (in case of a row-major
// matrix) or the row index (in case of a column-major matrix) must be a multiple of the number
// of values inside the intrinsic element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void StaticMatrix<Type,M,N,SO>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= NN , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   stream( &v_[i*NN+j], value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO2 >     // Storage order of the right-hand side dense matrix
inline typename DisableIf< typename StaticMatrix<Type,M,N,SO>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   StaticMatrix<Type,M,N,SO>::assign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         v_[i*NN+j] = (~rhs)(i,j);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO2 >     // Storage order of the right-hand side dense matrix
inline typename EnableIf< typename StaticMatrix<Type,M,N,SO>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   StaticMatrix<Type,M,N,SO>::assign( const DenseMatrix<MT,SO2>& rhs )
{
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; j+=IT::size ) {
         store( &v_[i*NN+j], (~rhs).load(i,j) );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,SO>::assign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t i=0UL; i<M; ++i )
      for( RhsConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i*NN+element->index()] = element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,SO>::assign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t j=0UL; j<N; ++j )
      for( RhsConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()*NN+j] = element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO2 >     // Storage order of the right-hand side dense matrix
inline typename DisableIf< typename StaticMatrix<Type,M,N,SO>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   StaticMatrix<Type,M,N,SO>::addAssign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         v_[i*NN+j] += (~rhs)(i,j);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO2 >     // Storage order of the right-hand side dense matrix
inline typename EnableIf< typename StaticMatrix<Type,M,N,SO>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   StaticMatrix<Type,M,N,SO>::addAssign( const DenseMatrix<MT,SO2>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; j+=IT::size ) {
         store( &v_[i*NN+j], load( &v_[i*NN+j] ) + (~rhs).load(i,j) );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,SO>::addAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t i=0UL; i<M; ++i )
      for( RhsConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i*NN+element->index()] += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,SO>::addAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t j=0UL; j<N; ++j )
      for( RhsConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()*NN+j] += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO2 >     // Storage order of the right-hand side dense matrix
inline typename DisableIf< typename StaticMatrix<Type,M,N,SO>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   StaticMatrix<Type,M,N,SO>::subAssign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         v_[i*NN+j] -= (~rhs)(i,j);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO2 >     // Storage order of the right-hand side dense matrix
inline typename EnableIf< typename StaticMatrix<Type,M,N,SO>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   StaticMatrix<Type,M,N,SO>::subAssign( const DenseMatrix<MT,SO2>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; j+=IT::size ) {
         store( &v_[i*NN+j], load( &v_[i*NN+j] ) - (~rhs).load(i,j) );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,SO>::subAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t i=0UL; i<M; ++i )
      for( RhsConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i*NN+element->index()] -= element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,SO>::subAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t j=0UL; j<N; ++j )
      for( RhsConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()*NN+j] -= element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of StaticMatrix for column-major matrices.
// \ingroup static_matrix
//
// This specialization of StaticMatrix adapts the class template to the requirements of
// column-major matrices.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
class StaticMatrix<Type,M,N,true> : public DenseMatrix< StaticMatrix<Type,M,N,true>, true >
{
 private:
   //**Type definitions****************************************************************************
   typedef IntrinsicTrait<Type>  IT;  //!< Intrinsic trait for the matrix element type.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Alignment adjustment
   enum { MM = M + ( IT::size - ( M % IT::size ) ) % IT::size };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef StaticMatrix<Type,M,N,true>   This;            //!< Type of this StaticMatrix instance.
   typedef This                          ResultType;      //!< Result type for expression template evaluations.
   typedef StaticMatrix<Type,M,N,false>  OppositeType;    //!< Result type with opposite storage order for expression template evaluations.
   typedef StaticMatrix<Type,N,M,false>  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef Type                          ElementType;     //!< Type of the matrix elements.
   typedef typename IT::Type             IntrinsicType;   //!< Intrinsic type of the matrix elements.
   typedef const Type&                   ReturnType;      //!< Return type for expression template evaluations.
   typedef const This&                   CompositeType;   //!< Data type for composite expression templates.
   typedef Type&                         Reference;       //!< Reference to a non-constant matrix value.
   typedef const Type&                   ConstReference;  //!< Reference to a constant matrix value.
   typedef DenseIterator<Type>           Iterator;        //!< Iterator over non-constant elements.
   typedef DenseIterator<const Type>     ConstIterator;   //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for intrinsic optimization.
   /*! The \a vectorizable compilation flag indicates whether expressions the matrix is involved
       in can be optimized via intrinsics. In case the element type of the matrix is a vectorizable
       data type, the \a vectorizable compilation flag is set to \a true, otherwise it is set to
       \a false. */
   enum { vectorizable = IsVectorizable<Type>::value };

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline StaticMatrix();
   explicit inline StaticMatrix( const Type& init );

   template< typename Other > explicit inline StaticMatrix( size_t m, size_t n, const Other* array );
   template< typename Other > explicit inline StaticMatrix( const Other (&array)[M][N] );

                                       inline StaticMatrix( const StaticMatrix& m );
   template< typename Other, bool SO > inline StaticMatrix( const StaticMatrix<Other,M,N,SO>&  m );
   template< typename MT   , bool SO > inline StaticMatrix( const Matrix<MT,SO>& m );

   inline StaticMatrix( const Type& v1, const Type& v2 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7, const Type& v8 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7, const Type& v8, const Type& v9 );
   inline StaticMatrix( const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5,
                        const Type& v6, const Type& v7, const Type& v8, const Type& v9, const Type& v10 );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Type*          data  ();
   inline const Type*    data  () const;
   inline Type*          data  ( size_t j );
   inline const Type*    data  ( size_t j ) const;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   template< typename Other >
   inline StaticMatrix& operator=( const Other (&array)[M][N] );

                                       inline StaticMatrix& operator= ( const Type& set );
                                       inline StaticMatrix& operator= ( const StaticMatrix& rhs );
   template< typename Other, bool SO > inline StaticMatrix& operator= ( const StaticMatrix<Other,M,N,SO>&  rhs );
   template< typename MT   , bool SO > inline StaticMatrix& operator= ( const Matrix<MT,SO>& rhs );
   template< typename MT   , bool SO > inline StaticMatrix& operator+=( const Matrix<MT,SO>& rhs );
   template< typename MT   , bool SO > inline StaticMatrix& operator-=( const Matrix<MT,SO>& rhs );
   template< typename MT   , bool SO > inline StaticMatrix& operator*=( const Matrix<MT,SO>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, StaticMatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, StaticMatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t        rows() const;
                              inline size_t        columns() const;
                              inline size_t        spacing() const;
                              inline size_t        capacity() const;
                              inline size_t        capacity( size_t j ) const;
                              inline size_t        nonZeros() const;
                              inline size_t        nonZeros( size_t j ) const;
                              inline void          reset();
                              inline void          reset( size_t i );
                              inline StaticMatrix& transpose();
   template< typename Other > inline StaticMatrix& scale( const Other& scalar );
                              inline void          swap( StaticMatrix& m ) /* throw() */;
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IsColumnMajorMatrix<MT>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IntrinsicTrait<Type>::addition &&
                     IsColumnMajorMatrix<MT>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IntrinsicTrait<Type>::subtraction &&
                     IsColumnMajorMatrix<MT>::value };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool isAligned() const;

   inline IntrinsicType load  ( size_t i, size_t j ) const;
   inline IntrinsicType loadu ( size_t i, size_t j ) const;
   inline void          store ( size_t i, size_t j, const IntrinsicType& value );
   inline void          storeu( size_t i, size_t j, const IntrinsicType& value );
   inline void          stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT, bool SO >
   inline typename DisableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT, bool SO >
   inline typename EnableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT > inline void assign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,false>& rhs );

   template< typename MT, bool SO >
   inline typename DisableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT, bool SO >
   inline typename EnableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT > inline void addAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,false>& rhs );

   template< typename MT, bool SO >
   inline typename DisableIf< VectorizedSubAssign<MT> >::Type
      subAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT, bool SO >
   inline typename EnableIf< VectorizedSubAssign<MT> >::Type
      subAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT > inline void subAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   AlignedArray<Type,MM*N> v_;  //!< The statically allocated matrix elements.
                                /*!< Access to the matrix elements is gained via the
                                     function call operator. */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_STATIC_ASSERT( MM % IT::size == 0UL );
   BLAZE_STATIC_ASSERT( MM >= M );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The default constructor for StaticMatrix.
//
// All matrix elements are initialized to the default value (i.e. 0 for integral data types).
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix()
{
   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<MM*N; ++i )
         v_[i] = Type();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a homogenous initialization of all elements.
//
// \param init Initial value for all matrix elements.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& init )
{
   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] = init;

      if( IsNumeric<Type>::value ) {
         for( size_t i=M; i<MM; ++i )
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array initialization of all matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param array Dynamic array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a dynamic array:

   \code
   using blaze::columnMajor;

   int* array = new int[6];
   // ... Initialization of the dynamic array
   blaze::StaticMatrix<int,3,4,columnMajor> v( array, 2UL, 3UL );
   delete[] array;
   \endcode

// The matrix is initialized with the values from the given array. Missing values are initialized
// with default values. In case the specified number of rows and/or columns exceeds the maximum
// number of rows/column of the static matrix (i.e. \m is larger than M or \a n is larger than N),
// a \a std::invalid_argument exception is thrown.\n
// Note that it is expected that the given \a array has at least \a m by \a n elements. Providing
// an array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the initialization array
inline StaticMatrix<Type,M,N,true>::StaticMatrix( size_t m, size_t n, const Other* array )
{
   if( m > M || n > N )
      throw std::invalid_argument( "Invalid setup of static matrix" );

   for( size_t j=0UL; j<n; ++j ) {
      for( size_t i=0UL; i<m; ++i )
         v_[i+j*MM] = array[i+j*m];

      if( IsNumeric<Type>::value ) {
         for( size_t i=m; i<MM; ++i )
            v_[i+j*MM] = Type();
      }
   }

   for( size_t j=n; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] = Type();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array initialization of all matrix elements.
//
// \param array \f$ M \times N \f$ dimensional array for the initialization.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a static array:

   \code
   using blaze::columnMajor;

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::StaticMatrix<int,3,3,columnMajor> A( init );
   \endcode

// The matrix is initialized with the values from the given array. Missing values are initialized
// with default values (as e.g. the value 6 in the example).
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the initialization array
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Other (&array)[M][N] )
{
   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] = array[i][j];

      if( IsNumeric<Type>::value ) {
         for( size_t i=M; i<MM; ++i )
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The copy constructor for StaticMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined in order to enable/facilitate NRV optimization.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const StaticMatrix& m )
{
   for( size_t i=0UL; i<MM*N; ++i )
      v_[i] = m.v_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion constructor from different StaticMatrix instances.
//
// \param m Matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , size_t M        // Number of rows
        , size_t N >      // Number of columns
template< typename Other  // Data type of the foreign matrix
        , bool SO >       // Storage order of the foreign matrix
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const StaticMatrix<Other,M,N,SO>& m )
{
   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] = m(i,j);

      if( IsNumeric<Type>::value ) {
         for( size_t i=M; i<MM; ++i )
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
// \exception std::invalid_argument Invalid setup of static matrix.
//
// This constructor initializes the static matrix from the given matrix. In case the size of
// the given matrix does not match the size of the static matrix (i.e. the number of rows is
// not M or the number of columns is not N), a \a std::invalid_argument exception is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the foreign matrix
        , bool SO >      // Storage order of the foreign matrix
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Matrix<MT,SO>& m )
{
   using blaze::assign;

   if( (~m).rows() != M || (~m).columns() != N )
      throw std::invalid_argument( "Invalid setup of static matrix" );

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=( IsSparseMatrix<MT>::value )?( 0UL ):( M ); i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }

   assign( *this, ~m );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 2 \f$ and \f$ 2 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 2 \f$
// and \f$ 2 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{2}{c}}
                     1 & 2 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,2,true> A( 1, 2 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2 )
{
   BLAZE_STATIC_ASSERT( M*N == 2UL );

   // Initialization of a 2x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
   }

   // Initialization of a 1x2 matrix
   else {
      v_[0UL] = v1;
      v_[ MM] = v2;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 3 \f$ and \f$ 3 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 3 \f$
// and \f$ 3 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{3}{c}}
                     1 & 2 & 3 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,3,true> A( 1, 2, 3 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3 )
{
   BLAZE_STATIC_ASSERT( M*N == 3UL );

   // Initialization of a 3x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
   }

   // Initialization of a 1x3 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 3 \f$, \f$ 2 \times 2 \f$, and \f$ 3 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 4 \f$,
// \f$ 2 \times 2 \f$, and \f$ 4 \times 1 \f$ matrix. The following example demonstrates this by
// creating the matrix

                     \f[\left(\begin{array}{*{2}{c}}
                     1 & 3 \\
                     2 & 4 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,2,true> A( 1, 2, 3, 4 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4 )
{
   BLAZE_STATIC_ASSERT( M*N == 4UL );

   // Initialization of a 4x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
   }

   // Initialization of a 2x2 matrix
   else if( N == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[MM    ] = v3;
      v_[MM+1UL] = v4;
   }

   // Initialization of a 1x4 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 5 \f$ and \f$ 5 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 5 \f$
// and \f$ 5 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{5}{c}}
                     1 & 2 & 3 & 4 & 5 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,5,true> A( 1, 2, 3, 4, 5 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4, const Type& v5 )
{
   BLAZE_STATIC_ASSERT( M*N == 5UL );

   // Initialization of a 5x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
   }

   // Initialization of a 1x5 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
      v_[4UL*MM] = v5;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 6 \f$, \f$ 2 \times 3 \f$, \f$ 3 \times 2 \f$, and
//        \f$ 6 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 6 \f$,
// \f$ 2 \times 3 \f$, \f$ 3 \times 2 \f$, and \f$ 6 \times 1 \f$ matrix. The following example
// demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{3}{c}}
                     1 & 3 & 5 \\
                     2 & 4 & 6 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,3,true> A( 1, 2, 3, 4, 5, 6 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4, const Type& v5, const Type& v6 )
{
   BLAZE_STATIC_ASSERT( M*N == 6UL );

   // Initialization of a 6x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
   }

   // Initialization of a 3x2 matrix
   else if( N == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[   2UL] = v3;
      v_[MM    ] = v4;
      v_[MM+1UL] = v5;
      v_[MM+2UL] = v6;
   }

   // Initialization of a 2x3 matrix
   else if( N == 3UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[    MM    ] = v3;
      v_[    MM+1UL] = v4;
      v_[2UL*MM    ] = v5;
      v_[2UL*MM+1UL] = v6;
   }

   // Initialization of a 1x6 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
      v_[4UL*MM] = v5;
      v_[5UL*MM] = v6;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 7 \f$ and \f$ 7 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 7 \f$
// and \f$ 7 \times 1 \f$ matrix. The following example demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{7}{c}}
                     1 & 2 & 3 & 4 & 5 & 6 & 7 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,1,7,true> A( 1, 2, 3, 4, 5, 6, 7 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4, const Type& v5, const Type& v6,
                                                  const Type& v7 )
{
   BLAZE_STATIC_ASSERT( M*N == 7UL );

   // Initialization of a 7x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
   }

   // Initialization of a 1x7 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
      v_[4UL*MM] = v5;
      v_[5UL*MM] = v6;
      v_[6UL*MM] = v7;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 8 \f$, \f$ 2 \times 4 \f$, \f$ 4 \times 2 \f$, and
//        \f$ 8 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
// \param v8 The initializer for the eigth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 8 \f$,
// \f$ 2 \times 4 \f$, \f$ 4 \times 2 \f$, and \f$ 8 \times 1 \f$ matrix. The following example
// demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{4}{c}}
                     1 & 3 & 5 & 7 \\
                     2 & 4 & 6 & 8 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,4,true> A( 1, 2, 3, 4, 5, 6, 7, 8 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4, const Type& v5, const Type& v6,
                                                  const Type& v7, const Type& v8 )
{
   BLAZE_STATIC_ASSERT( M*N == 8UL );

   // Initialization of a 8x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
      v_[7UL] = v8;
   }

   // Initialization of a 4x2 matrix
   else if( N == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[   2UL] = v3;
      v_[   3UL] = v4;
      v_[MM    ] = v5;
      v_[MM+1UL] = v6;
      v_[MM+2UL] = v7;
      v_[MM+3UL] = v8;
   }

   // Initialization of a 2x4 matrix
   else if( N == 4UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[    MM    ] = v3;
      v_[    MM+1UL] = v4;
      v_[2UL*MM    ] = v5;
      v_[2UL*MM+1UL] = v6;
      v_[3UL*MM    ] = v7;
      v_[3UL*MM+1UL] = v8;
   }

   // Initialization of a 1x8 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
      v_[4UL*MM] = v5;
      v_[5UL*MM] = v6;
      v_[6UL*MM] = v7;
      v_[7UL*MM] = v8;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 9 \f$, \f$ 3 \times 3 \f$, and \f$ 9 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
// \param v8 The initializer for the eigth matrix element.
// \param v9 The initializer for the ninth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 9 \f$,
// \f$ 3 \times 3 \f$, and \f$ 9 \times 1 \f$ matrix. The following example demonstrates this by
// creating the matrix

                     \f[\left(\begin{array}{*{3}{c}}
                     1 & 4 & 7 \\
                     2 & 5 & 8 \\
                     3 & 6 & 9 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,3,3,true> A( 1, 2, 3, 4, 5, 6, 7, 8, 9 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4, const Type& v5, const Type& v6,
                                                  const Type& v7, const Type& v8, const Type& v9 )
{
   BLAZE_STATIC_ASSERT( M*N == 9UL );

   // Initialization of a 9x1 matrix
   if( N == 1 ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
      v_[7UL] = v8;
      v_[8UL] = v9;
   }

   // Initialization of a 3x3 matrix
   else if( N == 3UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[       2UL] = v3;
      v_[    MM    ] = v4;
      v_[    MM+1UL] = v5;
      v_[    MM+2UL] = v6;
      v_[2UL*MM    ] = v7;
      v_[2UL*MM+1UL] = v8;
      v_[2UL*MM+2UL] = v9;
   }

   // Initialization of a 1x9 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
      v_[4UL*MM] = v5;
      v_[5UL*MM] = v6;
      v_[6UL*MM] = v7;
      v_[7UL*MM] = v8;
      v_[8UL*MM] = v9;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for \f$ 1 \times 10 \f$, \f$ 2 \times 5 \f$, \f$ 5 \times 2 \f$, and
// \f$ 10 \times 1 \f$ matrices.
//
// \param v1 The initializer for the first matrix element.
// \param v2 The initializer for the second matrix element.
// \param v3 The initializer for the third matrix element.
// \param v4 The initializer for the fourth matrix element.
// \param v5 The initializer for the fifth matrix element.
// \param v6 The initializer for the sixth matrix element.
// \param v7 The initializer for the seventh matrix element.
// \param v8 The initializer for the eigth matrix element.
// \param v9 The initializer for the ninth matrix element.
// \param v10 The initializer for the tenth matrix element.
//
// This constructor offers the option to directly initialize a newly created \f$ 1 \times 10 \f$,
// \f$ 2 \times 5 \f$, \f$ 5 \times 2 \f$, and \f$ 10 \times 1 \f$ matrix. The following example
// demonstrates this by creating the matrix

                     \f[\left(\begin{array}{*{5}{c}}
                     1 & 3 & 5 & 7 & 9  \\
                     2 & 4 & 6 & 8 & 10 \\
                     \end{array}\right)\f]:

   \code
   blaze::StaticMatrix<int,2,5,true> A( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
   \endcode
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>::StaticMatrix( const Type& v1, const Type& v2, const Type& v3,
                                                  const Type& v4, const Type& v5, const Type& v6,
                                                  const Type& v7, const Type& v8, const Type& v9,
                                                  const Type& v10 )
{
   BLAZE_STATIC_ASSERT( M*N == 10UL );

   // Initialization of a 10x1 matrix
   if( N == 1UL ) {
      v_[0UL] = v1;
      v_[1UL] = v2;
      v_[2UL] = v3;
      v_[3UL] = v4;
      v_[4UL] = v5;
      v_[5UL] = v6;
      v_[6UL] = v7;
      v_[7UL] = v8;
      v_[8UL] = v9;
      v_[9UL] = v10;
   }

   // Initialization of a 5x2 matrix
   else if( N == 2UL ) {
      v_[   0UL] = v1;
      v_[   1UL] = v2;
      v_[   2UL] = v3;
      v_[   3UL] = v4;
      v_[   4UL] = v5;
      v_[MM    ] = v6;
      v_[MM+1UL] = v7;
      v_[MM+2UL] = v8;
      v_[MM+3UL] = v9;
      v_[MM+4UL] = v10;
   }

   // Initialization of a 2x5 matrix
   else if( N == 5UL ) {
      v_[       0UL] = v1;
      v_[       1UL] = v2;
      v_[    MM    ] = v3;
      v_[    MM+1UL] = v4;
      v_[2UL*MM    ] = v5;
      v_[2UL*MM+1UL] = v6;
      v_[3UL*MM    ] = v7;
      v_[3UL*MM+1UL] = v8;
      v_[4UL*MM    ] = v9;
      v_[4UL*MM+1UL] = v10;
   }

   // Initialization of a 1x10 matrix
   else {
      v_[   0UL] = v1;
      v_[    MM] = v2;
      v_[2UL*MM] = v3;
      v_[3UL*MM] = v4;
      v_[4UL*MM] = v5;
      v_[5UL*MM] = v6;
      v_[6UL*MM] = v7;
      v_[7UL*MM] = v8;
      v_[8UL*MM] = v9;
      v_[9UL*MM] = v10;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=M; i<MM; ++i ) {
            v_[i+j*MM] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::Reference
   StaticMatrix<Type,M,N,true>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i<M, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<N, "Invalid column access index" );
   return v_[i+j*MM];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return Reference-to-const to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::ConstReference
   StaticMatrix<Type,M,N,true>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i<M, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<N, "Invalid column access index" );
   return v_[i+j*MM];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the static matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The static matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline Type* StaticMatrix<Type,M,N,true>::data()
{
   return v_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the static matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The static matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline const Type* StaticMatrix<Type,M,N,true>::data() const
{
   return v_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline Type* StaticMatrix<Type,M,N,true>::data( size_t j )
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return v_ + j*MM;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline const Type* StaticMatrix<Type,M,N,true>::data( size_t j ) const
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return v_ + j*MM;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of column \a j.
//
// \param j The column index.
// \return Iterator to the first element of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::Iterator
   StaticMatrix<Type,M,N,true>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return Iterator( v_ + j*MM );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of column \a j.
//
// \param j The column index.
// \return Iterator to the first element of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::ConstIterator
   StaticMatrix<Type,M,N,true>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*MM );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of column \a j.
//
// \param j The column index.
// \return Iterator to the first element of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::ConstIterator
   StaticMatrix<Type,M,N,true>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*MM );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last element of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::Iterator
   StaticMatrix<Type,M,N,true>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return Iterator( v_ + j*MM + M );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last element of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::ConstIterator
   StaticMatrix<Type,M,N,true>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*MM + M );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last element of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::ConstIterator
   StaticMatrix<Type,M,N,true>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < N, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*MM + M );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array assignment to all matrix elements.
//
// \param array \f$ M \times N \f$ dimensional array for the assignment.
// \return Reference to the assigned matrix.
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::rowMajor;

   const real init[3][3] = { { 1, 2, 3 },
                             { 4, 5 },
                             { 7, 8, 9 } };
   blaze::StaticMatrix<int,3UL,3UL,rowMajor> A;
   A = init;
   \endcode

// The matrix is assigned the values from the given array. Missing values are initialized with
// default values (as e.g. the value 6 in the example).
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the initialization array
inline StaticMatrix<Type,M,N,true>&
   StaticMatrix<Type,M,N,true>::operator=( const Other (&array)[M][N] )
{
   for( size_t j=0UL; j<N; ++j )
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] = array[i][j];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Homogenous assignment to all matrix elements.
//
// \param set Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>&
   StaticMatrix<Type,M,N,true>::operator=( const Type& set )
{
   for( size_t j=0UL; j<N; ++j )
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] = set;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for StaticMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// Explicit definition of a copy assignment operator for performance reasons.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>&
   StaticMatrix<Type,M,N,true>::operator=( const StaticMatrix& rhs )
{
   using blaze::assign;

   assign( *this, ~rhs );
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different StaticMatrix instances.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
*/
template< typename Type   // Data type of the matrix
        , size_t M        // Number of rows
        , size_t N >      // Number of columns
template< typename Other  // Data type of the foreign matrix
        , bool SO >       // Storage order of the foreign matrix
inline StaticMatrix<Type,M,N,true>&
   StaticMatrix<Type,M,N,true>::operator=( const StaticMatrix<Other,M,N,SO>& rhs )
{
   using blaze::assign;

   assign( *this, ~rhs );
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid assignment to static matrix.
//
// This constructor initializes the matrix as a copy of the given matrix. In case the
// number of rows of the given matrix is not M or the number of columns is not N, a
// \a std::invalid_argument exception is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,true>& StaticMatrix<Type,M,N,true>::operator=( const Matrix<MT,SO>& rhs )
{
   using blaze::assign;

   if( (~rhs).rows() != M || (~rhs).columns() != N )
      throw std::invalid_argument( "Invalid assignment to static matrix" );

   if( (~rhs).canAlias( this ) ) {
      StaticMatrix tmp( ~rhs );
      swap( tmp );
   }
   else {
      if( IsSparseMatrix<MT>::value )
         reset();
      assign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,true>& StaticMatrix<Type,M,N,true>::operator+=( const Matrix<MT,SO>& rhs )
{
   using blaze::addAssign;

   if( (~rhs).rows() != M || (~rhs).columns() != N )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      StaticMatrix tmp( ~rhs );
      addAssign( *this, tmp );
   }
   else {
      addAssign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,true>& StaticMatrix<Type,M,N,true>::operator-=( const Matrix<MT,SO>& rhs )
{
   using blaze::subAssign;

   if( (~rhs).rows() != M || (~rhs).columns() != N )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      StaticMatrix tmp( ~rhs );
      subAssign( *this, tmp );
   }
   else {
      subAssign( *this, ~rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline StaticMatrix<Type,M,N,true>& StaticMatrix<Type,M,N,true>::operator*=( const Matrix<MT,SO>& rhs )
{
   if( M != N || (~rhs).rows() != M || (~rhs).columns() != M )
      throw std::invalid_argument( "Matrix sizes do not match" );

   StaticMatrix tmp( *this * (~rhs) );
   return this->operator=( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication between a matrix and
//        a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, StaticMatrix<Type,M,N,true> >::Type&
   StaticMatrix<Type,M,N,true>::operator*=( Other rhs )
{
   using blaze::assign;

   assign( *this, (*this) * rhs );
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a matrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the matrix.
//
// \b Note: A division by zero is only checked by an user assert.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, StaticMatrix<Type,M,N,true> >::Type&
   StaticMatrix<Type,M,N,true>::operator/=( Other rhs )
{
   using blaze::assign;

   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   assign( *this, (*this) / rhs );
   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current number of rows of the matrix.
//
// \return The number of rows of the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::rows() const
{
   return M;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::columns() const
{
   return N;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two column, i.e. the total number
// of elements of a column.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::spacing() const
{
   return MM;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::capacity() const
{
   return MM*N;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::capacity( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return MM;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the total number of non-zero elements in the matrix
//
// \return The number of non-zero elements in the dense matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<N; ++j )
      for( size_t i=0UL; i<M; ++i )
         if( !isDefault( v_[i+j*MM] ) )
            ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline size_t StaticMatrix<Type,M,N,true>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( (j+1UL)*MM );
   size_t nonzeros( 0UL );

   for( size_t i=j*MM; i<iend; ++i )
      if( !isDefault( v_[i] ) )
         ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline void StaticMatrix<Type,M,N,true>::reset()
{
   using blaze::reset;

   for( size_t j=0UL; j<N; ++j )
      for( size_t i=0UL; i<M; ++i )
         reset( v_[i+j*MM] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function reset the values in the specified column to their default value. Note that
// the capacity of the column remains unchanged.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline void StaticMatrix<Type,M,N,true>::reset( size_t j )
{
   using blaze::reset;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   for( size_t i=0UL; i<M; ++i )
      reset( v_[i+j*MM] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Transposing the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline StaticMatrix<Type,M,N,true>& StaticMatrix<Type,M,N,true>::transpose()
{
   using std::swap;

   for( size_t j=1UL; j<N; ++j )
      for( size_t i=0UL; i<j; ++i )
         swap( v_[i+j*MM], v_[j+i*MM] );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A*=s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the matrix.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the scalar value
inline StaticMatrix<Type,M,N,true>&
   StaticMatrix<Type,M,N,true>::scale( const Other& scalar )
{
   for( size_t j=0UL; j<N; ++j )
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*MM] *= scalar;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Swapping the contents of two static matrices.
//
// \param m The matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline void StaticMatrix<Type,M,N,true>::swap( StaticMatrix& m ) /* throw() */
{
   using std::swap;

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         swap( v_[i+j*MM], m(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the matrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address can alias with the matrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the foreign expression
inline bool StaticMatrix<Type,M,N,true>::canAlias( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the matrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address is aliased with the matrix. In contrast
// to the conAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
template< typename Other >  // Data type of the foreign expression
inline bool StaticMatrix<Type,M,N,true>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the matrix is properly aligned in memory.
//
// \return \a true in case the matrix is aligned, \a false if not.
//
// This function returns whether the matrix is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of each column of the matrix are guaranteed to conform to
// the alignment restrictions of the element type \a Type.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline bool StaticMatrix<Type,M,N,true>::isAligned() const
{
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an aligned load of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::IntrinsicType
   StaticMatrix<Type,M,N,true>::load( size_t i, size_t j ) const
{
   using blaze::load;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= MM , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );

   return load( &v_[i+j*MM] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded intrinsic element.
//
// This function performs an unaligned load of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline typename StaticMatrix<Type,M,N,true>::IntrinsicType
   StaticMatrix<Type,M,N,true>::loadu( size_t i, size_t j ) const
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= MM , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );

   return loadu( &v_[i+j*MM] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned store of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline void StaticMatrix<Type,M,N,true>::store( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::store;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= MM , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );

   store( &v_[i+j*MM], value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an unaligned store of a specific intrinsic element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the intrinsic element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline void StaticMatrix<Type,M,N,true>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= MM, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N , "Invalid column access index" );

   storeu( &v_[i+j*MM], value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of an intrinsic element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The intrinsic element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific intrinsic element of the
// dense matrix. // The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a multiple
// of the number of values inside the intrinsic element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
inline void StaticMatrix<Type,M,N,true>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  M  , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= MM , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  N  , "Invalid column access index" );

   stream( &v_[i+j*MM], value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO >      // Storage order of the right-hand side dense matrix
inline typename DisableIf< typename StaticMatrix<Type,M,N,true>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   StaticMatrix<Type,M,N,true>::assign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         v_[i+j*MM] = (~rhs)(i,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO >      // Storage order of the right-hand side dense matrix
inline typename EnableIf< typename StaticMatrix<Type,M,N,true>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   StaticMatrix<Type,M,N,true>::assign( const DenseMatrix<MT,SO>& rhs )
{
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; i+=IT::size ) {
         store( &v_[i+j*MM], (~rhs).load(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,true>::assign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t j=0UL; j<N; ++j )
      for( RhsConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()+j*MM] = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,true>::assign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t i=0UL; i<M; ++i )
      for( RhsConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i+element->index()*MM] = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO >      // Storage order of the right-hand side dense matrix
inline typename DisableIf< typename StaticMatrix<Type,M,N,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   StaticMatrix<Type,M,N,true>::addAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         v_[i+j*MM] += (~rhs)(i,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO >      // Storage order of the right-hand side dense matrix
inline typename EnableIf< typename StaticMatrix<Type,M,N,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   StaticMatrix<Type,M,N,true>::addAssign( const DenseMatrix<MT,SO>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; i+=IT::size ) {
         store( &v_[i+j*MM], load( &v_[i+j*MM] ) + (~rhs).load(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,true>::addAssign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t j=0UL; j<N; ++j )
      for( RhsConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()+j*MM] += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,true>::addAssign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t i=0UL; i<M; ++i )
      for( RhsConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i+element->index()*MM] += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO >      // Storage order of the right-hand side dense matrix
inline typename DisableIf< typename StaticMatrix<Type,M,N,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   StaticMatrix<Type,M,N,true>::subAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         v_[i+j*MM] -= (~rhs)(i,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT    // Type of the right-hand side dense matrix
        , bool SO >      // Storage order of the right-hand side dense matrix
inline typename EnableIf< typename StaticMatrix<Type,M,N,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   StaticMatrix<Type,M,N,true>::subAssign( const DenseMatrix<MT,SO>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; i+=IT::size ) {
         store( &v_[i+j*MM], load( &v_[i+j*MM] ) - (~rhs).load(i,j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,true>::subAssign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t j=0UL; j<N; ++j )
      for( RhsConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()+j*MM] -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
template< typename MT >  // Type of the right-hand side sparse matrix
inline void StaticMatrix<Type,M,N,true>::subAssign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() == M && (~rhs).columns() == N, "Invalid matrix size" );

   typedef typename MT::ConstIterator  RhsConstIterator;

   for( size_t i=0UL; i<M; ++i )
      for( RhsConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i+element->index()*MM] -= element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  UNDEFINED CLASS TEMPLATE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of StaticMatrix for zero columns.
// \ingroup static_matrix
//
// This specialization of the StaticMatrix class template is left undefined and therefore
// prevents the instantiation for zero columns.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , bool SO >      // Storage order
class StaticMatrix<Type,M,0UL,SO>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of StaticMatrix for zero rows.
// \ingroup static_matrix
//
// This specialization of the StaticMatrix class template is left undefined and therefore
// prevents the instantiation for zero rows.
*/
template< typename Type  // Data type of the matrix
        , size_t N       // Number of columns
        , bool SO >      // Storage order
class StaticMatrix<Type,0UL,N,SO>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of StaticMatrix for 0 rows and 0 columns.
// \ingroup static_matrix
//
// This specialization of the StaticMatrix class template is left undefined and therefore
// prevents the instantiation for 0 rows and 0 columns.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
class StaticMatrix<Type,0UL,0UL,SO>;
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  STATICMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name StaticMatrix operators */
//@{
template< typename Type, size_t M, size_t N, bool SO >
inline void reset( StaticMatrix<Type,M,N,SO>& m );

template< typename Type, size_t M, size_t N, bool SO >
inline void clear( StaticMatrix<Type,M,N,SO>& m );

template< typename Type, size_t M, size_t N, bool SO >
inline bool isDefault( const StaticMatrix<Type,M,N,SO>& m );

template< typename Type, size_t M, size_t N, bool SO >
inline void swap( StaticMatrix<Type,M,N,SO>& a, StaticMatrix<Type,M,N,SO>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given static matrix.
// \ingroup static_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void reset( StaticMatrix<Type,M,N,SO>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given static matrix.
// \ingroup static_matrix
//
// \param m The matrix to be cleared.
// \return void
//
// Clearing a static matrix is equivalent to resetting it via the reset() function.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void clear( StaticMatrix<Type,M,N,SO>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given static matrix is in default state.
// \ingroup static_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline bool isDefault( const StaticMatrix<Type,M,N,SO>& m )
{
   if( SO == rowMajor ) {
      for( size_t i=0UL; i<M; ++i )
         for( size_t j=0UL; j<N; ++j )
            if( !isDefault( m(i,j) ) ) return false;
   }
   else {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=0UL; i<M; ++i )
            if( !isDefault( m(i,j) ) ) return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two static matrices.
// \ingroup static_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline void swap( StaticMatrix<Type,M,N,SO>& a, StaticMatrix<Type,M,N,SO>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct AddTrait< StaticMatrix<T1,M,N,SO>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticMatrix< typename AddTrait<T1,T2>::Type, M, N, SO >  Type;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct AddTrait< StaticMatrix<T1,M,N,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   typedef StaticMatrix< typename AddTrait<T1,T2>::Type, M, N, false >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct SubTrait< StaticMatrix<T1,M,N,SO>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticMatrix< typename SubTrait<T1,T2>::Type, M, N, SO >  Type;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct SubTrait< StaticMatrix<T1,M,N,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   typedef StaticMatrix< typename SubTrait<T1,T2>::Type, M, N, false >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct MultTrait< StaticMatrix<T1,M,N,SO>, T2 >
{
   typedef StaticMatrix< typename MultTrait<T1,T2>::Type, M, N, SO >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};

template< typename T1, typename T2, size_t M, size_t N, bool SO >
struct MultTrait< T1, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticMatrix< typename MultTrait<T1,T2>::Type, M, N, SO >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T1 );
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct MultTrait< StaticMatrix<T1,M,N,SO>, StaticVector<T2,N,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, M, false >  Type;
};

template< typename T1, size_t M, typename T2, size_t N, bool SO >
struct MultTrait< StaticVector<T1,M,true>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, N, true >  Type;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2, size_t L >
struct MultTrait< StaticMatrix<T1,M,N,SO>, HybridVector<T2,L,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, M, false >  Type;
};

template< typename T1, size_t L, typename T2, size_t M, size_t N, bool SO >
struct MultTrait< HybridVector<T1,L,true>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, N, true >  Type;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct MultTrait< StaticMatrix<T1,M,N,SO>, DynamicVector<T2,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, M, false >  Type;
};

template< typename T1, typename T2, size_t M, size_t N, bool SO >
struct MultTrait< DynamicVector<T1,true>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, N, true >  Type;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct MultTrait< StaticMatrix<T1,M,N,SO>, CompressedVector<T2,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, M, false >  Type;
};

template< typename T1, typename T2, size_t M, size_t N, bool SO >
struct MultTrait< CompressedVector<T1,true>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, N, true >  Type;
};

template< typename T1, size_t M, size_t K, bool SO1, typename T2, size_t N, bool SO2 >
struct MultTrait< StaticMatrix<T1,M,K,SO1>, StaticMatrix<T2,K,N,SO2> >
{
   typedef StaticMatrix< typename MultTrait<T1,T2>::Type, M, N, SO1 >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct DivTrait< StaticMatrix<T1,M,N,SO>, T2 >
{
   typedef StaticMatrix< typename DivTrait<T1,T2>::Type, M, N, SO >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MATHTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct MathTrait< StaticMatrix<T1,M,N,SO>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticMatrix< typename MathTrait<T1,T2>::HighType, M, N, SO >  HighType;
   typedef StaticMatrix< typename MathTrait<T1,T2>::LowType , M, N, SO >  LowType;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO >
struct SubmatrixTrait< StaticMatrix<T1,M,N,SO> >
{
   typedef DynamicMatrix<T1,SO>  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO >
struct RowTrait< StaticMatrix<T1,M,N,SO> >
{
   typedef StaticVector<T1,N,true>  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, size_t M, size_t N, bool SO >
struct ColumnTrait< StaticMatrix<T1,M,N,SO> >
{
   typedef StaticVector<T1,M,false>  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
