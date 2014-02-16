//=================================================================================================
/*!
//  \file blaze/math/dense/DynamicMatrix.h
//  \brief Header file for the implementation of a dynamic MxN matrix
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

#ifndef _BLAZE_MATH_DENSE_DYNAMICMATRIX_H_
#define _BLAZE_MATH_DENSE_DYNAMICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <stdexcept>
#include <blaze/math/dense/DenseIterator.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/Functions.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/smp/DenseMatrix.h>
#include <blaze/math/smp/SparseMatrix.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MathTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Restrict.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/system/Streaming.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Null.h>
#include <blaze/util/Template.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dynamic_matrix DynamicMatrix
// \ingroup dense_matrix
*/
/*!\brief Efficient implementation of a dynamic \f$ M \times N \f$ matrix.
// \ingroup dynamic_matrix
//
// The DynamicMatrix class template is the representation of an arbitrary sized matrix with
// \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. The type of the elements
// and the storage order of the matrix can be specified via the two template parameters:

   \code
   template< typename Type, bool SO >
   class DynamicMatrix;
   \endcode

//  - Type: specifies the type of the matrix elements. DynamicMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
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

// The use of DynamicMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combination of row-major and
// column-major dense and sparse matrices with fitting element types. The following example gives
// an impression of the use of DynamicMatrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   DynamicMatrix<double,rowMajor> A( 2, 3 );  // Default constructed, non-initialized, row-major 2x3 matrix
   A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;  // Initialization of the first row
   A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;  // Initialization of the second row

   DynamicMatrix<float,columnMajor> B( 2, 3 );  // Default constructed column-major single precision 2x3 matrix
   B(0,0) = 1.0; B(0,1) = 3.0; B(0,2) = 5.0;    // Initialization of the first row
   B(1,0) = 2.0; B(1,1) = 4.0; B(1,2) = 6.0;    // Initialization of the second row

   CompressedMatrix<float> C( 2, 3 );        // Empty row-major sparse single precision matrix
   DynamicMatrix<float>    D( 3, 2, 4.0F );  // Directly, homogeneously initialized single precision 3x2 matrix

   DynamicMatrix<double,rowMajor>    E( A );  // Creation of a new row-major matrix as a copy of A
   DynamicMatrix<double,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;     // Matrix addition and assignment to a row-major matrix
   E = A - C;     // Matrix subtraction and assignment to a column-major matrix
   F = A * D;     // Matrix multiplication between two matrices of different element types

   A *= 2.0;      // In-place scaling of matrix A
   E  = 2.0 * B;  // Scaling of matrix B
   F  = D * 2.0;  // Scaling of matrix D

   E += A - B;    // Addition assignment
   E -= A + C;    // Subtraction assignment
   F *= A * D;    // Multiplication assignment
   \endcode
*/
template< typename Type                    // Data type of the matrix
        , bool SO = defaultStorageOrder >  // Storage order
class DynamicMatrix : public DenseMatrix< DynamicMatrix<Type,SO>, SO >
{
 private:
   //**Type definitions****************************************************************************
   typedef IntrinsicTrait<Type>  IT;  //!< Intrinsic trait for the matrix element type.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DynamicMatrix<Type,SO>     This;            //!< Type of this DynamicMatrix instance.
   typedef This                       ResultType;      //!< Result type for expression template evaluations.
   typedef DynamicMatrix<Type,!SO>    OppositeType;    //!< Result type with opposite storage order for expression template evaluations.
   typedef DynamicMatrix<Type,!SO>    TransposeType;   //!< Transpose type for expression template evaluations.
   typedef Type                       ElementType;     //!< Type of the matrix elements.
   typedef typename IT::Type          IntrinsicType;   //!< Intrinsic type of the matrix elements.
   typedef const Type&                ReturnType;      //!< Return type for expression template evaluations.
   typedef const This&                CompositeType;   //!< Data type for composite expression templates.
   typedef Type&                      Reference;       //!< Reference to a non-constant matrix value.
   typedef const Type&                ConstReference;  //!< Reference to a constant matrix value.
   typedef DenseIterator<Type>        Iterator;        //!< Iterator over non-constant elements.
   typedef DenseIterator<const Type>  ConstIterator;   //!< Iterator over constant elements.
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
   enum { smpAssignable = 1 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
                              explicit inline DynamicMatrix();
                              explicit inline DynamicMatrix( size_t m, size_t n );
                              explicit inline DynamicMatrix( size_t m, size_t n, const Type& init );
   template< typename Other > explicit inline DynamicMatrix( size_t m, size_t n, const Other* array );

   template< typename Other, size_t M, size_t N >
   explicit inline DynamicMatrix( const Other (&array)[M][N] );

                                     inline DynamicMatrix( const DynamicMatrix& m );
   template< typename MT, bool SO2 > inline DynamicMatrix( const Matrix<MT,SO2>& m );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~DynamicMatrix();
   //@}
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
   template< typename Other, size_t M, size_t N >
   inline DynamicMatrix& operator=( const Other (&array)[M][N] );

                                     inline DynamicMatrix& operator= ( Type set );
                                     inline DynamicMatrix& operator= ( const DynamicMatrix& set );
   template< typename MT, bool SO2 > inline DynamicMatrix& operator= ( const Matrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline DynamicMatrix& operator+=( const Matrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline DynamicMatrix& operator-=( const Matrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline DynamicMatrix& operator*=( const Matrix<MT,SO2>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DynamicMatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DynamicMatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t         rows() const;
                              inline size_t         columns() const;
                              inline size_t         spacing() const;
                              inline size_t         capacity() const;
                              inline size_t         capacity( size_t i ) const;
                              inline size_t         nonZeros() const;
                              inline size_t         nonZeros( size_t i ) const;
                              inline void           reset();
                              inline void           reset( size_t i );
                              inline void           clear();
                                     void           resize ( size_t m, size_t n, bool preserve=true );
                              inline void           extend ( size_t m, size_t n, bool preserve=true );
                              inline void           reserve( size_t elements );
                              inline DynamicMatrix& transpose();
   template< typename Other > inline DynamicMatrix& scale( Other scalar );
                              inline void           swap( DynamicMatrix& m ) /* throw() */;
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value };
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
                     IntrinsicTrait<Type>::addition };
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
                     IntrinsicTrait<Type>::subtraction };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load  ( size_t i, size_t j ) const;
   inline IntrinsicType loadu ( size_t i, size_t j ) const;
   inline void          store ( size_t i, size_t j, const IntrinsicType& value );
   inline void          storeu( size_t i, size_t j, const IntrinsicType& value );
   inline void          stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT >
   inline typename DisableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT >
   inline typename EnableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT > inline void assign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT >
   inline typename DisableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT >
   inline typename EnableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT > inline void addAssign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT >
   inline typename DisableIf< VectorizedSubAssign<MT> >::Type
      subAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT >
   inline typename EnableIf< VectorizedSubAssign<MT> >::Type
      subAssign( const DenseMatrix<MT,SO>& rhs );

   template< typename MT > inline void subAssign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,!SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t adjustColumns( size_t minColumns ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;                //!< The current number of rows of the matrix.
   size_t n_;                //!< The current number of columns of the matrix.
   size_t nn_;               //!< The alignment adjusted number of columns.
   size_t capacity_;         //!< The maximum capacity of the matrix.
   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated matrix elements.
                             /*!< Access to the matrix elements is gained via the subscript or
                                  function call operator. In case of row-major order the memory
                                  layout of the elements is
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
/*!\brief The default constructor for DynamicMatrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>::DynamicMatrix()
   : m_       ( 0UL  )  // The current number of rows of the matrix
   , n_       ( 0UL  )  // The current number of columns of the matrix
   , nn_      ( 0UL  )  // The alignment adjusted number of columns
   , capacity_( 0UL  )  // The maximum capacity of the matrix
   , v_       ( NULL )  // The matrix elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$. No element initialization is performed!
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
//
// \b Note: This constructor is only responsible to allocate the required dynamic memory. No
//          element initialization is performed!
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>::DynamicMatrix( size_t m, size_t n )
   : m_       ( m )                            // The current number of rows of the matrix
   , n_       ( n )                            // The current number of columns of the matrix
   , nn_      ( adjustColumns( n ) )           // The alignment adjusted number of columns
   , capacity_( m_*nn_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<m_; ++i ) {
         for( size_t j=n_; j<nn_; ++j )
            v_[i*nn_+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param init The initial value of the matrix elements.
//
// All matrix elements are initialized with the specified value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>::DynamicMatrix( size_t m, size_t n, const Type& init )
   : m_       ( m )                            // The current number of rows of the matrix
   , n_       ( n )                            // The current number of columns of the matrix
   , nn_      ( adjustColumns( n ) )           // The alignment adjusted number of columns
   , capacity_( m_*nn_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   for( size_t i=0UL; i<m; ++i ) {
      for( size_t j=0UL; j<n_; ++j )
         v_[i*nn_+j] = init;

      if( IsNumeric<Type>::value ) {
         for( size_t j=n_; j<nn_; ++j )
            v_[i*nn_+j] = Type();
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

   int* array = new int[20];
   // ... Initialization of the dynamic array
   blaze::DynamicMatrix<int,rowMajor> v( 4UL, 5UL, array );
   delete[] array;
   \endcode

// The matrix is sized accoring to the given size of the array and initialized with the values
// from the given array. Note that it is expected that the given \a array has at least \a m by
// \a n elements. Providing an array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the initialization array
inline DynamicMatrix<Type,SO>::DynamicMatrix( size_t m, size_t n, const Other* array )
   : m_       ( m )                            // The current number of rows of the matrix
   , n_       ( n )                            // The current number of columns of the matrix
   , nn_      ( adjustColumns( n ) )           // The alignment adjusted number of columns
   , capacity_( m_*nn_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   for( size_t i=0UL; i<m; ++i ) {
      for( size_t j=0UL; j<n; ++j )
         v_[i*nn_+j] = array[i*n+j];

      if( IsNumeric<Type>::value ) {
         for( size_t j=n; j<nn_; ++j )
            v_[i*nn_+j] = Type();
      }
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
   blaze::DynamicMatrix<int,rowMajor> A( init );
   \endcode

// The matrix is sized accoring to the size of the array and initialized with the values from
// the given array. Missing values are initialized with default values (as e.g. the value 6 in
// the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO >       // Storage order
template< typename Other  // Data type of the initialization array
        , size_t M        // Number of rows of the initialization array
        , size_t N >      // Number of columns of the initialization array
inline DynamicMatrix<Type,SO>::DynamicMatrix( const Other (&array)[M][N] )
   : m_       ( M )                            // The current number of rows of the matrix
   , n_       ( N )                            // The current number of columns of the matrix
   , nn_      ( adjustColumns( N ) )           // The alignment adjusted number of columns
   , capacity_( m_*nn_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j )
         v_[i*nn_+j] = array[i][j];

      if( IsNumeric<Type>::value ) {
         for( size_t j=N; j<nn_; ++j )
            v_[i*nn_+j] = Type();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for DynamicMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>::DynamicMatrix( const DynamicMatrix& m )
   : m_       ( m.m_  )                        // The current number of rows of the matrix
   , n_       ( m.n_  )                        // The current number of columns of the matrix
   , nn_      ( m.nn_ )                        // The alignment adjusted number of columns
   , capacity_( m_*nn_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   BLAZE_INTERNAL_ASSERT( capacity_ <= m.capacity_, "Invalid capacity estimation" );

   for( size_t i=0UL; i<capacity_; ++i )
      v_[i] = m.v_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
template< typename MT    // Type of the foreign matrix
        , bool SO2 >     // Storage order of the foreign matrix
inline DynamicMatrix<Type,SO>::DynamicMatrix( const Matrix<MT,SO2>& m )
   : m_       ( (~m).rows() )                  // The current number of rows of the matrix
   , n_       ( (~m).columns() )               // The current number of columns of the matrix
   , nn_      ( adjustColumns( n_ ) )          // The alignment adjusted number of columns
   , capacity_( m_*nn_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<m_; ++i ) {
         for( size_t j=( IsSparseMatrix<MT>::value )?( 0UL ):( n_ ); j<nn_; ++j )
            v_[i*nn_+j] = Type();
      }
   }

   smpAssign( *this, ~m );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for DynamicMatrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>::~DynamicMatrix()
{
   deallocate( v_ );
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
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::Reference
   DynamicMatrix<Type,SO>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i*nn_+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::ConstReference
   DynamicMatrix<Type,SO>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i*nn_+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Type* DynamicMatrix<Type,SO>::data()
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const Type* DynamicMatrix<Type,SO>::data() const
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
        , bool SO >      // Storage order
inline Type* DynamicMatrix<Type,SO>::data( size_t i )
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return v_ + i*nn_;
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
        , bool SO >      // Storage order
inline const Type* DynamicMatrix<Type,SO>::data( size_t i ) const
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return v_ + i*nn_;
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::Iterator
   DynamicMatrix<Type,SO>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return Iterator( v_ + i*nn_ );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::ConstIterator
   DynamicMatrix<Type,SO>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::ConstIterator
   DynamicMatrix<Type,SO>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::Iterator
   DynamicMatrix<Type,SO>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return Iterator( v_ + i*nn_ + n_ );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::ConstIterator
   DynamicMatrix<Type,SO>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ + n_ );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::ConstIterator
   DynamicMatrix<Type,SO>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ + n_ );
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
   blaze::DynamicMatrix<int,rowMajor> A;
   A = init;
   \endcode

// The matrix is resized accoring to the size of the array and assigned the values of the given
// array. Missing values are initialized with default values (as e.g. the value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO >       // Storage order
template< typename Other  // Data type of the initialization array
        , size_t M        // Number of rows of the initialization array
        , size_t N >      // Number of columns of the initialization array
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator=( const Other (&array)[M][N] )
{
   resize( M, N, false );

   for( size_t i=0UL; i<M; ++i )
      for( size_t j=0UL; j<N; ++j )
         v_[i*nn_+j] = array[i][j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Homogenous assignment to all matrix elements.
//
// \param rhs Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator=( Type rhs )
{
   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         v_[i*nn_+j] = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DynamicMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator=( const DynamicMatrix& rhs )
{
   if( &rhs == this ) return *this;

   resize( rhs.m_, rhs.n_, false );
   smpAssign( *this, ~rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator=( const Matrix<MT,SO2>& rhs )
{
   if( (~rhs).canAlias( this ) ) {
      DynamicMatrix tmp( ~rhs );
      swap( tmp );
   }
   else {
      resize( (~rhs).rows(), (~rhs).columns(), false );
      if( IsSparseMatrix<MT>::value )
         reset();
      smpAssign( *this, ~rhs );
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
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator+=( const Matrix<MT,SO2>& rhs )
{
   if( (~rhs).rows() != m_ || (~rhs).columns() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      typename MT::ResultType tmp( ~rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, ~rhs );
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
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator-=( const Matrix<MT,SO2>& rhs )
{
   if( (~rhs).rows() != m_ || (~rhs).columns() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      typename MT::ResultType tmp( ~rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, ~rhs );
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
        , bool SO >      // Storage order
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::operator*=( const Matrix<MT,SO2>& rhs )
{
   if( (~rhs).rows() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   DynamicMatrix tmp( *this * (~rhs) );
   swap( tmp );

   return *this;
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
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DynamicMatrix<Type,SO> >::Type&
   DynamicMatrix<Type,SO>::operator*=( Other rhs )
{
   smpAssign( *this, (*this) * rhs );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a matrix by a scalar value
//        (\f$ A/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the matrix.
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DynamicMatrix<Type,SO> >::Type&
   DynamicMatrix<Type,SO>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   smpAssign( *this, (*this) / rhs );
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
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::rows() const
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::columns() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows/columns.
//
// \return The spacing between the beginning of two rows/columns.
//
// This function returns the spacing between the beginning of two rows/columns, i.e. the
// total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::spacing() const
{
   return nn_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::capacity() const
{
   return capacity_;
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
template< typename Type  // Data type of the sparse matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::capacity( size_t i ) const
{
   UNUSED_PARAMETER( i );
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return nn_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the total number of non-zero elements in the matrix
//
// \return The number of non-zero elements in the dense matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         if( !isDefault( v_[i*nn_+j] ) )
            ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( (i+1UL)*nn_ );
   size_t nonzeros( 0UL );

   for( size_t j=i*nn_; j<jend; ++j )
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
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::reset()
{
   using blaze::reset;

   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         reset( v_[i*nn_+j] );
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
template< typename Type  // Data type of the sparse matrix
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::reset( size_t i )
{
   using blaze::reset;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   for( size_t j=0UL; j<n_; ++j )
      reset( v_[i*nn_+j] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the \f$ M \times N \f$ matrix.
//
// \return void
//
// After the clear() function, the size of the matrix is 0.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::clear()
{
   resize( 0UL, 0UL, false );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the matrix.
//
// \param m The new number of rows of the matrix.
// \param n The new number of columns of the matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Therefore this function invalidates all existing views (rows, columns, ...) on
// the matrix. Additionally, the resize operation potentially changes all matrix elements.
// In order to preserve the old matrix values, the \a preserve flag can be set to \a true.
// However, new matrix elements are not initialized!\n
// The following example illustrates the resize operation of a \f$ 2 \times 4 \f$ matrix to a
// \f$ 4 \times 2 \f$ matrix. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & 3 & 4 \\
                              5 & 6 & 7 & 8 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              5 & 6 \\
                              x & x \\
                              x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
void DynamicMatrix<Type,SO>::resize( size_t m, size_t n, bool preserve )
{
   using blaze::min;

   if( m == m_ && n == n_ ) return;

   const size_t nn( adjustColumns( n ) );

   if( preserve )
   {
      Type* BLAZE_RESTRICT v = allocate<Type>( m*nn );
      const size_t min_m( min( m, m_ ) );
      const size_t min_n( min( n, n_ ) );

      for( size_t i=0UL; i<min_m; ++i )
         for( size_t j=0UL; j<min_n; ++j )
            v[i*nn+j] = v_[i*nn_+j];

      std::swap( v_, v );
      deallocate( v );
      capacity_ = m*nn;
   }
   else if( m*nn > capacity_ ) {
      Type* BLAZE_RESTRICT v = allocate<Type>( m*nn );
      std::swap( v_, v );
      deallocate( v );
      capacity_ = m*nn;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t i=0UL; i<m; ++i )
         for( size_t j=n; j<nn; ++j )
            v_[i*nn+j] = Type();
   }

   m_  = m;
   n_  = n;
   nn_ = nn;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the matrix.
//
// \param m Number of additional rows.
// \param n Number of additional columns.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function increases the matrix size by \a m rows and \a n columns. During this operation,
// new dynamic memory may be allocated in case the capacity of the matrix is too small. Therefore
// this function potentially changes all matrix elements. In order to preserve the old matrix
// values, the \a preserve flag can be set to \a true. However, new matrix elements are not
// initialized!
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::extend( size_t m, size_t n, bool preserve )
{
   resize( m_+m, n_+n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the matrix.
//
// \param elements The new minimum capacity of the sparse matrix.
// \return void
//
// This function increases the capacity of the sparse matrix to at least \a elements elements.
// The current values of the matrix elements are preserved.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::reserve( size_t elements )
{
   if( elements > capacity_ )
   {
      // Allocating a new array
      Type* BLAZE_RESTRICT tmp = allocate<Type>( elements );

      // Initializing the new array
      std::copy( v_, v_+capacity_, tmp );

      if( IsNumeric<Type>::value ) {
         for( size_t i=capacity_; i<elements; ++i )
            tmp[i] = Type();
      }

      // Replacing the old array
      std::swap( tmp, v_ );
      deallocate( tmp );
      capacity_ = elements;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transposing the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::transpose()
{
   DynamicMatrix tmp( trans(*this) );
   swap( tmp );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the matrix.
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the scalar value
inline DynamicMatrix<Type,SO>& DynamicMatrix<Type,SO>::scale( Other scalar )
{
   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         v_[i*nn_+j] *= scalar;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
//
// \param m The matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::swap( DynamicMatrix& m ) /* throw() */
{
   std::swap( m_ , m.m_  );
   std::swap( n_ , m.n_  );
   std::swap( nn_, m.nn_ );
   std::swap( capacity_, m.capacity_ );
   std::swap( v_ , m.v_  );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Adjusting the number columns of the matrix according to its data type \a Type.
//
// \param minColumns The minimum necessary number of columns.
// \return The adjusted number of columns.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t DynamicMatrix<Type,SO>::adjustColumns( size_t minColumns ) const
{
   if( IsNumeric<Type>::value )
      return minColumns + ( IT::size - ( minColumns % IT::size ) ) % IT::size;
   else return minColumns;
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
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,SO>::canAlias( const Other* alias ) const
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
// This function returns whether the given address is aliased with the vector. In contrast
// to the conAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , bool SO >         // Storage order
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,SO>::isAliased( const Other* alias ) const
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
        , bool SO >      // Storage order
inline bool DynamicMatrix<Type,SO>::isAligned() const
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix can be used in SMP assignments.
//
// \return \a true in case the matrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the matrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline bool DynamicMatrix<Type,SO>::canSMPAssign() const
{
   return ( rows() > OPENMP_DMATASSIGN_THRESHOLD );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::IntrinsicType
   DynamicMatrix<Type,SO>::load( size_t i, size_t j ) const
{
   using blaze::load;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= nn_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   return load( v_+i*nn_+j );
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
        , bool SO >      // Storage order
inline typename DynamicMatrix<Type,SO>::IntrinsicType
   DynamicMatrix<Type,SO>::loadu( size_t i, size_t j ) const
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= nn_, "Invalid column access index" );

   return loadu( v_+i*nn_+j );
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
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::store( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::store;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= nn_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   store( v_+i*nn_+j, value );
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
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= nn_, "Invalid column access index" );

   storeu( v_+i*nn_+j, value );
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
        , bool SO >      // Storage order
inline void DynamicMatrix<Type,SO>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + IT::size <= nn_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % IT::size == 0UL, "Invalid column access index" );

   stream( v_+i*nn_+j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a row-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DynamicMatrix<Type,SO>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   DynamicMatrix<Type,SO>::assign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jend( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jend, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jend; j+=2UL ) {
         v_[i*nn_+j    ] = (~rhs)(i,j    );
         v_[i*nn_+j+1UL] = (~rhs)(i,j+1UL);
      }
      if( jend < n_ ) {
         v_[i*nn_+jend] = (~rhs)(i,jend);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the assignment of a row-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DynamicMatrix<Type,SO>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   DynamicMatrix<Type,SO>::assign( const DenseMatrix<MT,SO>& rhs )
{
   using blaze::store;
   using blaze::stream;

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   if( useStreaming && m_*n_ > ( cacheSize / ( sizeof(Type) * 3UL ) ) && !(~rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<m_; ++i )
         for( size_t j=0UL; j<n_; j+=IT::size )
            stream( v_+i*nn_+j, (~rhs).load(i,j) );
   }
   else
   {
      const size_t jend( n_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

      for( size_t i=0UL; i<m_; ++i ) {
         typename MT::ConstIterator it( (~rhs).begin(i) );
         for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
            store( v_+i*nn_+j             , it.load() ); it += IT::size;
            store( v_+i*nn_+j+IT::size    , it.load() ); it += IT::size;
            store( v_+i*nn_+j+IT::size*2UL, it.load() ); it += IT::size;
            store( v_+i*nn_+j+IT::size*3UL, it.load() ); it += IT::size;
         }
         for( size_t j=jend; j<n_; j+=IT::size, it+=IT::size ) {
            store( v_+i*nn_+j, it.load() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a column-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO>::assign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               v_[i*nn_+j] = (~rhs)(i,j);
            }
         }
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO>::assign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i*nn_+element->index()] = element->value();
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO>::assign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()*nn_+j] = element->value();
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DynamicMatrix<Type,SO>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   DynamicMatrix<Type,SO>::addAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jend( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jend, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jend; j+=2UL ) {
         v_[i*nn_+j    ] += (~rhs)(i,j    );
         v_[i*nn_+j+1UL] += (~rhs)(i,j+1UL);
      }
      if( jend < n_ ) {
         v_[i*nn_+jend] += (~rhs)(i,jend);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the addition assignment of a row-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DynamicMatrix<Type,SO>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   DynamicMatrix<Type,SO>::addAssign( const DenseMatrix<MT,SO>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   const size_t jend( n_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      typename MT::ConstIterator it( (~rhs).begin(i) );
      for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
         store( v_+i*nn_+j             , load( v_+i*nn_+j              ) + it.load() ); it += IT::size;
         store( v_+i*nn_+j+IT::size    , load( v_+i*nn_+j+IT::size     ) + it.load() ); it += IT::size;
         store( v_+i*nn_+j+IT::size*2UL, load( v_+i*nn_+j+IT::size*2UL ) + it.load() ); it += IT::size;
         store( v_+i*nn_+j+IT::size*3UL, load( v_+i*nn_+j+IT::size*3UL ) + it.load() ); it += IT::size;
      }
      for( size_t j=jend; j<n_; j+=IT::size, it+=IT::size ) {
         store( v_+i*nn_+j, load( v_+i*nn_+j ) + it.load() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO>::addAssign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               v_[i*nn_+j] += (~rhs)(i,j);
            }
         }
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO>::addAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i*nn_+element->index()] += element->value();
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO>::addAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()*nn_+j] += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline typename DisableIf< typename DynamicMatrix<Type,SO>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   DynamicMatrix<Type,SO>::subAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t jend( n_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % 2UL ) ) == jend, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jend; j+=2UL ) {
         v_[i*nn_+j    ] -= (~rhs)(i,j    );
         v_[i*nn_+j+1UL] -= (~rhs)(i,j+1UL);
      }
      if( jend < n_ ) {
         v_[i*nn_+jend] -= (~rhs)(i,jend);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a row-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline typename EnableIf< typename DynamicMatrix<Type,SO>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   DynamicMatrix<Type,SO>::subAssign( const DenseMatrix<MT,SO>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   const size_t jend( n_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( n_ - ( n_ % (IT::size*4UL) ) ) == jend, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      typename MT::ConstIterator it( (~rhs).begin(i) );
      for( size_t j=0UL; j<jend; j+=IT::size*4UL ) {
         store( v_+i*nn_+j             , load( v_+i*nn_+j              ) - it.load() ); it += IT::size;
         store( v_+i*nn_+j+IT::size    , load( v_+i*nn_+j+IT::size     ) - it.load() ); it += IT::size;
         store( v_+i*nn_+j+IT::size*2UL, load( v_+i*nn_+j+IT::size*2UL ) - it.load() ); it += IT::size;
         store( v_+i*nn_+j+IT::size*3UL, load( v_+i*nn_+j+IT::size*3UL ) - it.load() ); it += IT::size;
      }
      for( size_t j=jend; j<n_; j+=IT::size, it+=IT::size ) {
         store( v_+i*nn_+j, load( v_+i*nn_+j ) - it.load() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO>::subAssign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               v_[i*nn_+j] -= (~rhs)(i,j);
            }
         }
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO>::subAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( typename MT::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i*nn_+element->index()] -= element->value();
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
        , bool SO >      // Storage order
template< typename MT >  // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO>::subAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( typename MT::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()*nn_+j] -= element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of DynamicMatrix for column-major matrices.
// \ingroup dynamic_matrix
//
// This specialization of DynamicMatrix adapts the class template to the requirements of
// column-major matrices.
*/
template< typename Type >  // Data type of the matrix
class DynamicMatrix<Type,true> : public DenseMatrix< DynamicMatrix<Type,true>, true >
{
 private:
   //**Type definitions****************************************************************************
   typedef IntrinsicTrait<Type>  IT;  //!< Intrinsic trait for the matrix element type.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef DynamicMatrix<Type,true>   This;            //!< Type of this DynamicMatrix instance.
   typedef This                       ResultType;      //!< Result type for expression template evaluations.
   typedef DynamicMatrix<Type,false>  OppositeType;    //!< Result type with opposite storage order for expression template evaluations.
   typedef DynamicMatrix<Type,false>  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef Type                       ElementType;     //!< Type of the matrix elements.
   typedef typename IT::Type          IntrinsicType;   //!< Intrinsic type of the matrix elements.
   typedef const Type&                ReturnType;      //!< Return type for expression template evaluations.
   typedef const This&                CompositeType;   //!< Data type for composite expression templates.
   typedef Type&                      Reference;       //!< Reference to a non-constant matrix value.
   typedef const Type&                ConstReference;  //!< Reference to a constant matrix value.
   typedef DenseIterator<Type>        Iterator;        //!< Iterator over non-constant elements.
   typedef DenseIterator<const Type>  ConstIterator;   //!< Iterator over constant elements.
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
   enum { smpAssignable = 1 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
                              explicit inline DynamicMatrix();
                              explicit inline DynamicMatrix( size_t m, size_t n );
                              explicit inline DynamicMatrix( size_t m, size_t n, const Type& init );
   template< typename Other > explicit inline DynamicMatrix( size_t m, size_t n, const Other* array );

   template< typename Other, size_t M, size_t N >
   explicit inline DynamicMatrix( const Other (&array)[M][N] );

                                    inline DynamicMatrix( const DynamicMatrix& m );
   template< typename MT, bool SO > inline DynamicMatrix( const Matrix<MT,SO>& m );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~DynamicMatrix();
   //@}
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
   template< typename Other, size_t M, size_t N >
   inline DynamicMatrix& operator=( const Other (&array)[M][N] );

                                    inline DynamicMatrix& operator= ( Type set );
                                    inline DynamicMatrix& operator= ( const DynamicMatrix& set );
   template< typename MT, bool SO > inline DynamicMatrix& operator= ( const Matrix<MT,SO>& rhs );
   template< typename MT, bool SO > inline DynamicMatrix& operator+=( const Matrix<MT,SO>& rhs );
   template< typename MT, bool SO > inline DynamicMatrix& operator-=( const Matrix<MT,SO>& rhs );
   template< typename MT, bool SO > inline DynamicMatrix& operator*=( const Matrix<MT,SO>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DynamicMatrix >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, DynamicMatrix >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t         rows() const;
                              inline size_t         columns() const;
                              inline size_t         spacing() const;
                              inline size_t         capacity() const;
                              inline size_t         capacity( size_t j ) const;
                              inline size_t         nonZeros() const;
                              inline size_t         nonZeros( size_t j ) const;
                              inline void           reset();
                              inline void           reset( size_t j );
                              inline void           clear();
                                     void           resize ( size_t m, size_t n, bool preserve=true );
                              inline void           extend ( size_t m, size_t n, bool preserve=true );
                              inline void           reserve( size_t elements );
                              inline DynamicMatrix& transpose();
   template< typename Other > inline DynamicMatrix& scale( Other scalar );
                              inline void           swap( DynamicMatrix& m ) /* throw() */;
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedAddAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IntrinsicTrait<Type>::addition };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename MT >
   struct VectorizedSubAssign {
      enum { value = vectorizable && MT::vectorizable &&
                     IsSame<Type,typename MT::ElementType>::value &&
                     IntrinsicTrait<Type>::subtraction };
   };
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   inline bool isAligned   () const;
   inline bool canSMPAssign() const;

   inline IntrinsicType load  ( size_t i, size_t j ) const;
   inline IntrinsicType loadu ( size_t i, size_t j ) const;
   inline void          store ( size_t i, size_t j, const IntrinsicType& value );
   inline void          storeu( size_t i, size_t j, const IntrinsicType& value );
   inline void          stream( size_t i, size_t j, const IntrinsicType& value );

   template< typename MT >
   inline typename DisableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,true>& rhs );

   template< typename MT >
   inline typename EnableIf< VectorizedAssign<MT> >::Type
      assign( const DenseMatrix<MT,true>& rhs );

   template< typename MT > inline void assign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,false>& rhs );

   template< typename MT >
   inline typename DisableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,true>& rhs );

   template< typename MT >
   inline typename EnableIf< VectorizedAddAssign<MT> >::Type
      addAssign( const DenseMatrix<MT,true>& rhs );

   template< typename MT > inline void addAssign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,false>& rhs );

   template< typename MT >
   inline typename DisableIf< VectorizedSubAssign<MT> >::Type
      subAssign ( const DenseMatrix<MT,true>& rhs );

   template< typename MT >
   inline typename EnableIf< VectorizedSubAssign<MT> >::Type
      subAssign ( const DenseMatrix<MT,true>& rhs );

   template< typename MT > inline void subAssign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t adjustRows( size_t minRows ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;                //!< The current number of rows of the matrix.
   size_t mm_;               //!< The alignment adjusted number of rows.
   size_t n_;                //!< The current number of columns of the matrix.
   size_t capacity_;         //!< The maximum capacity of the matrix.
   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated matrix elements.
                             /*!< Access to the matrix elements is gained via the subscript or
                                  function call operator. In case of row-major order the memory
                                  layout of the elements is
                                  \f[\left(\begin{array}{*{5}{c}}
                                  0            & 1             & 2             & \cdots & N-1         \\
                                  N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                                  \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                                  M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
                                  \end{array}\right)\f]. */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
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
/*!\brief The default constructor for DynamicMatrix.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>::DynamicMatrix()
   : m_       ( 0UL  )  // The current number of rows of the matrix
   , mm_      ( 0UL  )  // The alignment adjusted number of rows
   , n_       ( 0UL  )  // The current number of columns of the matrix
   , capacity_( 0UL  )  // The maximum capacity of the matrix
   , v_       ( NULL )  // The matrix elements
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ m \times n \f$. No element initialization is performed!
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
//
// \b Note: This constructor is only responsible to allocate the required dynamic memory. No
//          element initialization is performed!
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>::DynamicMatrix( size_t m, size_t n )
   : m_       ( m )                            // The current number of rows of the matrix
   , mm_      ( adjustRows( m ) )              // The alignment adjusted number of rows
   , n_       ( n )                            // The current number of columns of the matrix
   , capacity_( mm_*n_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=m_; i<mm_; ++i ) {
            v_[i+j*mm_] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param init The initial value of the matrix elements.
//
// All matrix elements are initialized with the specified value.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>::DynamicMatrix( size_t m, size_t n, const Type& init )
   : m_       ( m )                            // The current number of rows of the matrix
   , mm_      ( adjustRows( m ) )              // The alignment adjusted number of rows
   , n_       ( n )                            // The current number of columns of the matrix
   , capacity_( mm_*n_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<m_; ++i )
         v_[i+j*mm_] = init;

      if( IsNumeric<Type>::value ) {
         for( size_t i=m_; i<mm_; ++i )
            v_[i+j*mm_] = Type();
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

   int* array = new int[20];
   // ... Initialization of the dynamic array
   blaze::DynamicMatrix<int,columnMajor> v( array, 5UL, 4UL );
   delete[] array;
   \endcode

// The matrix is sized accoring to the given size of the array and initialized with the values
// from the given array. Note that it is expected that the given \a array has at least \a m by
// \a n elements. Providing an array with less elements results in undefined behavior!
*/
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the initialization array
inline DynamicMatrix<Type,true>::DynamicMatrix( size_t m, size_t n, const Other* array )
   : m_       ( m )                            // The current number of rows of the matrix
   , mm_      ( adjustRows( m ) )              // The alignment adjusted number of rows
   , n_       ( n )                            // The current number of columns of the matrix
   , capacity_( mm_*n_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   for( size_t j=0UL; j<n; ++j ) {
      for( size_t i=0UL; i<m; ++i )
         v_[i+j*mm_] = array[i+j*m];

      if( IsNumeric<Type>::value ) {
         for( size_t i=m; i<mm_; ++i )
            v_[i+j*mm_] = Type();
      }
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
   blaze::DynamicMatrix<int,columnMajor> A( init );
   \endcode

// The matrix is sized according to the size of the array and initialized with the values from
// the given array. Missing values are initialized with default values (as e.g. the value 6 in
// the example).
*/
template< typename Type >  // Data type of the matrix
template< typename Other   // Data type of the initialization array
        , size_t M         // Number of rows of the initialization array
        , size_t N >       // Number of columns of the initialization array
inline DynamicMatrix<Type,true>::DynamicMatrix( const Other (&array)[M][N] )
   : m_       ( M )                            // The current number of rows of the matrix
   , mm_      ( adjustRows( M ) )              // The alignment adjusted number of rows
   , n_       ( N )                            // The current number of columns of the matrix
   , capacity_( mm_*n_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*mm_] = array[i][j];

      if( IsNumeric<Type>::value ) {
         for( size_t i=M; i<mm_; ++i )
            v_[i+j*mm_] = Type();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The copy constructor for DynamicMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>::DynamicMatrix( const DynamicMatrix& m )
   : m_       ( m.m_  )                        // The current number of rows of the matrix
   , mm_      ( m.mm_ )                        // The alignment adjusted number of rows
   , n_       ( m.n_  )                        // The current number of columns of the matrix
   , capacity_( mm_*n_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   BLAZE_INTERNAL_ASSERT( capacity_ <= m.capacity_, "Invalid capacity estimation" );

   for( size_t i=0UL; i<capacity_; ++i )
      v_[i] = m.v_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
*/
template< typename Type >  // Data type of the matrix
template< typename MT      // Type of the foreign matrix
        , bool SO >        // Storage order of the foreign matrix
inline DynamicMatrix<Type,true>::DynamicMatrix( const Matrix<MT,SO>& m )
   : m_       ( (~m).rows() )                  // The current number of rows of the matrix
   , mm_      ( adjustRows( m_ ) )             // The alignment adjusted number of rows
   , n_       ( (~m).columns() )               // The current number of columns of the matrix
   , capacity_( mm_*n_ )                       // The maximum capacity of the matrix
   , v_       ( allocate<Type>( capacity_ ) )  // The matrix elements
{
   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=( IsSparseMatrix<MT>::value )?( 0UL ):( m_ ); i<mm_; ++i ) {
            v_[i+j*mm_] = Type();
      }
   }

   smpAssign( *this, ~m );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The destructor for DynamicMatrix.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>::~DynamicMatrix()
{
   deallocate( v_ );
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
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::Reference
   DynamicMatrix<Type,true>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i+j*mm_];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::ConstReference
   DynamicMatrix<Type,true>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i+j*mm_];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type >  // Data type of the matrix
inline Type* DynamicMatrix<Type,true>::data()
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
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data.
*/
template< typename Type >  // Data type of the matrix
inline const Type* DynamicMatrix<Type,true>::data() const
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
template< typename Type >  // Data type of the matrix
inline Type* DynamicMatrix<Type,true>::data( size_t j )
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return v_ + j*mm_;
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
template< typename Type >  // Data type of the matrix
inline const Type* DynamicMatrix<Type,true>::data( size_t j ) const
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return v_ + j*mm_;
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::Iterator
   DynamicMatrix<Type,true>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return Iterator( v_ + j*mm_ );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::ConstIterator
   DynamicMatrix<Type,true>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::ConstIterator
   DynamicMatrix<Type,true>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::Iterator
   DynamicMatrix<Type,true>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return Iterator( v_ + j*mm_ + m_ );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::ConstIterator
   DynamicMatrix<Type,true>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ + m_ );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::ConstIterator
   DynamicMatrix<Type,true>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ + m_ );
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

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::DynamicMatrix<int,rowMajor> A;
   A = init;
   \endcode

// The matrix is resized accoring to the size of the array and assigned the values of the given
// array. Missing values are initialized with default values (as e.g. the value 6 in the example).
*/
template< typename Type >  // Data type of the matrix
template< typename Other   // Data type of the initialization array
        , size_t M         // Number of rows of the initialization array
        , size_t N >       // Number of columns of the initialization array
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator=( const Other (&array)[M][N] )
{
   resize( M, N, false );

   for( size_t j=0UL; j<N; ++j )
      for( size_t i=0UL; i<M; ++i )
         v_[i+j*mm_] = array[i][j];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Homogenous assignment to all matrix elements.
//
// \param rhs Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator=( Type rhs )
{
   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         v_[i+j*mm_] = rhs;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DynamicMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator=( const DynamicMatrix& rhs )
{
   if( &rhs == this ) return *this;

   resize( rhs.m_, rhs.n_, false );
   smpAssign( *this, ~rhs );

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
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type >  // Data type of the matrix
template< typename MT      // Type of the right-hand side matrix
        , bool SO >        // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator=( const Matrix<MT,SO>& rhs )
{
   if( (~rhs).canAlias( this ) ) {
      DynamicMatrix tmp( ~rhs );
      swap( tmp );
   }
   else {
      resize( (~rhs).rows(), (~rhs).columns(), false );
      if( IsSparseMatrix<MT>::value )
         reset();
      smpAssign( *this, ~rhs );
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
template< typename Type >  // Data type of the matrix
template< typename MT      // Type of the right-hand side matrix
        , bool SO >        // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator+=( const Matrix<MT,SO>& rhs )
{
   if( (~rhs).rows() != m_ || (~rhs).columns() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      typename MT::ResultType tmp( ~rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, ~rhs );
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
template< typename Type >  // Data type of the matrix
template< typename MT      // Type of the right-hand side matrix
        , bool SO >        // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator-=( const Matrix<MT,SO>& rhs )
{
   if( (~rhs).rows() != m_ || (~rhs).columns() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   if( (~rhs).canAlias( this ) ) {
      typename MT::ResultType tmp( ~rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, ~rhs );
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
template< typename Type >  // Data type of the matrix
template< typename MT      // Type of the right-hand side matrix
        , bool SO >        // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::operator*=( const Matrix<MT,SO>& rhs )
{
   if( (~rhs).rows() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   DynamicMatrix tmp( *this * (~rhs) );
   swap( tmp );

   return *this;
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
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DynamicMatrix<Type,true> >::Type&
   DynamicMatrix<Type,true>::operator*=( Other rhs )
{
   smpAssign( *this, (*this) * rhs );
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
*/
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, DynamicMatrix<Type,true> >::Type&
   DynamicMatrix<Type,true>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   smpAssign( *this, (*this) / rhs );
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
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::rows() const
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::columns() const
{
   return n_;
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
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::spacing() const
{
   return mm_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::capacity() const
{
   return capacity_;
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
template< typename Type >  // Data type of the sparse matrix
inline size_t DynamicMatrix<Type,true>::capacity( size_t j ) const
{
   UNUSED_PARAMETER( j );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return mm_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the total number of non-zero elements in the matrix
//
// \return The number of non-zero elements in the dense matrix.
*/
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         if( !isDefault( v_[i+j*mm_] ) )
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
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( (j+1UL)*mm_ );
   size_t nonzeros( 0UL );

   for( size_t i=j*mm_; i<iend; ++i )
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
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::reset()
{
   using blaze::reset;

   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         reset( v_[i+j*mm_] );
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
template< typename Type >  // Data type of the sparse matrix
inline void DynamicMatrix<Type,true>::reset( size_t j )
{
   using blaze::reset;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   for( size_t i=0UL; i<m_; ++i )
      reset( v_[i+j*mm_] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the \f$ M \times N \f$ matrix.
//
// \return void
//
// After the clear() function, the size of the matrix is 0.
*/
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::clear()
{
   resize( 0UL, 0UL, false );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Changing the size of the matrix.
//
// \param m The new number of rows of the matrix.
// \param n The new number of columns of the matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Therefore this function invalidates all existing views (rows, columns, ...) on
// the matrix. Additionally, the resize operation potentially changes all matrix elements.
// In order to preserve the old matrix values, the \a preserve flag can be set to \a true.
// However, new matrix elements are not initialized!\n
// The following example illustrates the resize operation of a \f$ 2 \times 4 \f$ matrix to a
// \f$ 4 \times 2 \f$ matrix. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & 3 & 4 \\
                              5 & 6 & 7 & 8 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              5 & 6 \\
                              x & x \\
                              x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type >  // Data type of the matrix
void DynamicMatrix<Type,true>::resize( size_t m, size_t n, bool preserve )
{
   using blaze::min;

   if( m == m_ && n == n_ ) return;

   const size_t mm( adjustRows( m ) );

   if( preserve )
   {
      Type* BLAZE_RESTRICT v = allocate<Type>( mm*n );
      const size_t min_m( min( m, m_ ) );
      const size_t min_n( min( n, n_ ) );

      for( size_t j=0UL; j<min_n; ++j )
         for( size_t i=0UL; i<min_m; ++i )
            v[i+j*mm] = v_[i+j*mm_];

      std::swap( v_, v );
      deallocate( v );
      capacity_ = mm*n;
   }
   else if( mm*n > capacity_ ) {
      Type* BLAZE_RESTRICT v = allocate<Type>( mm*n );
      std::swap( v_, v );
      deallocate( v );
      capacity_ = mm*n;
   }

   if( IsNumeric<Type>::value ) {
      for( size_t j=0UL; j<n; ++j )
         for( size_t i=m; i<mm; ++i )
            v_[i+j*mm] = Type();
   }

   m_  = m;
   mm_ = mm;
   n_  = n;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Extending the size of the matrix.
//
// \param m Number of additional rows.
// \param n Number of additional columns.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function increases the matrix size by \a m rows and \a n columns. During this operation,
// new dynamic memory may be allocated in case the capacity of the matrix is too small. Therefore
// this function potentially changes all matrix elements. In order to preserve the old matrix
// values, the \a preserve flag can be set to \a true. However, new matrix elements are not
// initialized!
*/
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::extend( size_t m, size_t n, bool preserve )
{
   resize( m_+m, n_+n, preserve );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the matrix.
//
// \param elements The new minimum capacity of the sparse matrix.
// \return void
//
// This function increases the capacity of the sparse matrix to at least \a elements elements.
// The current values of the matrix elements are preserved.
*/
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::reserve( size_t elements )
{
   if( elements > capacity_ )
   {
      // Allocating a new array
      Type* BLAZE_RESTRICT tmp = allocate<Type>( elements );

      // Initializing the new array
      std::copy( v_, v_+capacity_, tmp );

      if( IsNumeric<Type>::value ) {
         for( size_t i=capacity_; i<elements; ++i )
            tmp[i] = Type();
      }

      // Replacing the old array
      std::swap( tmp, v_ );
      deallocate( tmp );
      capacity_ = elements;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Transposing the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type >  // Data type of the matrix
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::transpose()
{
   DynamicMatrix tmp( trans(*this) );
   swap( tmp );
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the matrix.
*/
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the scalar value
inline DynamicMatrix<Type,true>& DynamicMatrix<Type,true>::scale( Other scalar )
{
   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         v_[i+j*mm_] *= scalar;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Swapping the contents of two matrices.
//
// \param m The matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::swap( DynamicMatrix& m ) /* throw() */
{
   std::swap( m_ , m.m_  );
   std::swap( mm_, m.mm_ );
   std::swap( n_ , m.n_  );
   std::swap( capacity_, m.capacity_ );
   std::swap( v_ , m.v_  );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Adjusting the number rows of the matrix according to its data type \a Type.
//
// \param minRows The minimum necessary number of rows.
// \return The adjusted number of rows.
*/
template< typename Type >  // Data type of the matrix
inline size_t DynamicMatrix<Type,true>::adjustRows( size_t minRows ) const
{
   if( IsNumeric<Type>::value )
      return minRows + ( IT::size - ( minRows % IT::size ) ) % IT::size;
   else return minRows;
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
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,true>::canAlias( const Other* alias ) const
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
// This function returns whether the given address is aliased with the vector. In contrast
// to the conAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,true>::isAliased( const Other* alias ) const
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
template< typename Type >  // Data type of the matrix
inline bool DynamicMatrix<Type,true>::isAligned() const
{
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the matrix can be used in SMP assignments.
//
// \return \a true in case the matrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the matrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename Type >  // Data type of the matrix
inline bool DynamicMatrix<Type,true>::canSMPAssign() const
{
   return ( columns() > OPENMP_DMATASSIGN_THRESHOLD );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::IntrinsicType
   DynamicMatrix<Type,true>::load( size_t i, size_t j ) const
{
   using blaze::load;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= mm_, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );

   return load( v_+i+j*mm_ );
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
template< typename Type >  // Data type of the matrix
inline typename DynamicMatrix<Type,true>::IntrinsicType
   DynamicMatrix<Type,true>::loadu( size_t i, size_t j ) const
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= mm_, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );

   return loadu( v_+i+j*mm_ );
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
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::store( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::store;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= mm_, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );

   store( v_+i+j*mm_, value );
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
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::storeu( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= mm_, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );

   storeu( v_+i+j*mm_, value );
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
template< typename Type >  // Data type of the matrix
inline void DynamicMatrix<Type,true>::stream( size_t i, size_t j, const IntrinsicType& value )
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i            <  m_ , "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i + IT::size <= mm_, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( i % IT::size == 0UL, "Invalid row access index"    );
   BLAZE_INTERNAL_ASSERT( j            <  n_ , "Invalid column access index" );

   stream( v_+i+j*mm_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline typename DisableIf< typename DynamicMatrix<Type,true>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   DynamicMatrix<Type,true>::assign( const DenseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t iend( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == iend, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<iend; i+=2UL ) {
         v_[i    +j*mm_] = (~rhs)(i    ,j);
         v_[i+1UL+j*mm_] = (~rhs)(i+1UL,j);
      }
      if( iend < m_ ) {
         v_[iend+j*mm_] = (~rhs)(iend,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline typename EnableIf< typename DynamicMatrix<Type,true>::BLAZE_TEMPLATE VectorizedAssign<MT> >::Type
   DynamicMatrix<Type,true>::assign( const DenseMatrix<MT,true>& rhs )
{
   using blaze::store;
   using blaze::stream;

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   if( useStreaming && m_*n_ > ( cacheSize / ( sizeof(Type) * 3UL ) ) && !(~rhs).isAliased( this ) )
   {
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=0UL; i<m_; i+=IT::size )
            stream( v_+i+j*mm_, (~rhs).load(i,j) );
   }
   else
   {
      const size_t iend( m_ & size_t(-IT::size*4) );
      BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

      for( size_t j=0UL; j<n_; ++j ) {
         typename MT::ConstIterator it( (~rhs).begin(j) );
         for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
            store( v_+i+j*mm_             , it.load() ); it += IT::size;
            store( v_+i+j*mm_+IT::size    , it.load() ); it += IT::size;
            store( v_+i+j*mm_+IT::size*2UL, it.load() ); it += IT::size;
            store( v_+i+j*mm_+IT::size*3UL, it.load() ); it += IT::size;
         }
         for( size_t i=iend; i<m_; i+=IT::size, it+=IT::size ) {
            store( v_+i+j*mm_, it.load() );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true>::assign( const DenseMatrix<MT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               v_[i+j*mm_] = (~rhs)(i,j);
            }
         }
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true>::assign( const SparseMatrix<MT,true>& rhs )
{
   for( size_t j=0UL; j<(~rhs).columns(); ++j )
      for( typename MT::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()+j*mm_] = element->value();
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true>::assign( const SparseMatrix<MT,false>& rhs )
{
   for( size_t i=0UL; i<(~rhs).rows(); ++i )
      for( typename MT::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i+element->index()*mm_] = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline typename DisableIf< typename DynamicMatrix<Type,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   DynamicMatrix<Type,true>::addAssign( const DenseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t iend( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == iend, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<iend; i+=2UL ) {
         v_[i    +j*mm_] += (~rhs)(i    ,j);
         v_[i+1UL+j*mm_] += (~rhs)(i+1UL,j);
      }
      if( iend < m_ ) {
         v_[iend+j*mm_] += (~rhs)(iend,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline typename EnableIf< typename DynamicMatrix<Type,true>::BLAZE_TEMPLATE VectorizedAddAssign<MT> >::Type
   DynamicMatrix<Type,true>::addAssign( const DenseMatrix<MT,true>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   const size_t iend( m_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      typename MT::ConstIterator it( (~rhs).begin(j) );
      for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
         store( v_+i+j*mm_             , load( v_+i+j*mm_              ) + it.load() ); it += IT::size;
         store( v_+i+j*mm_+IT::size    , load( v_+i+j*mm_+IT::size     ) + it.load() ); it += IT::size;
         store( v_+i+j*mm_+IT::size*2UL, load( v_+i+j*mm_+IT::size*2UL ) + it.load() ); it += IT::size;
         store( v_+i+j*mm_+IT::size*3UL, load( v_+i+j*mm_+IT::size*3UL ) + it.load() ); it += IT::size;
      }
      for( size_t i=iend; i<m_; i+=IT::size, it+=IT::size ) {
         store( v_+i+j*mm_, load( v_+i+j*mm_ ) + it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true>::addAssign( const DenseMatrix<MT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               v_[i+j*mm_] += (~rhs)(i,j);
            }
         }
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true>::addAssign( const SparseMatrix<MT,true>& rhs )
{
   for( size_t j=0UL; j<(~rhs).columns(); ++j )
      for( typename MT::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()+j*mm_] += element->value();
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true>::addAssign( const SparseMatrix<MT,false>& rhs )
{
   for( size_t i=0UL; i<(~rhs).rows(); ++i )
      for( typename MT::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i+element->index()*mm_] += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline typename DisableIf< typename DynamicMatrix<Type,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   DynamicMatrix<Type,true>::subAssign( const DenseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   const size_t iend( m_ & size_t(-2) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % 2UL ) ) == iend, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<iend; i+=2UL ) {
         v_[i  +j*mm_] -= (~rhs)(i  ,j);
         v_[i+1+j*mm_] -= (~rhs)(i+1,j);
      }
      if( iend < m_ ) {
         v_[iend+j*mm_] -= (~rhs)(iend,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Intrinsic optimized implementation of the subtraction assignment of a column-major
//        dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline typename EnableIf< typename DynamicMatrix<Type,true>::BLAZE_TEMPLATE VectorizedSubAssign<MT> >::Type
   DynamicMatrix<Type,true>::subAssign( const DenseMatrix<MT,true>& rhs )
{
   using blaze::load;
   using blaze::store;

   BLAZE_INTERNAL_ASSERT( m_ == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (~rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   const size_t iend( m_ & size_t(-IT::size*4) );
   BLAZE_INTERNAL_ASSERT( ( m_ - ( m_ % (IT::size*4UL) ) ) == iend, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      typename MT::ConstIterator it( (~rhs).begin(j) );
      for( size_t i=0UL; i<iend; i+=IT::size*4UL ) {
         store( v_+i+j*mm_             , load( v_+i+j*mm_              ) - it.load() ); it += IT::size;
         store( v_+i+j*mm_+IT::size    , load( v_+i+j*mm_+IT::size     ) - it.load() ); it += IT::size;
         store( v_+i+j*mm_+IT::size*2UL, load( v_+i+j*mm_+IT::size*2UL ) - it.load() ); it += IT::size;
         store( v_+i+j*mm_+IT::size*3UL, load( v_+i+j*mm_+IT::size*3UL ) - it.load() ); it += IT::size;
      }
      for( size_t i=iend; i<m_; i+=IT::size, it+=IT::size ) {
         store( v_+i+j*mm_, load( v_+i+j*mm_ ) - it.load() );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true>::subAssign( const DenseMatrix<MT,false>& rhs )
{
   const size_t block( 16UL );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( ( n_<(jj+block) )?( n_ ):( jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( ( m_<(ii+block) )?( m_ ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               v_[i+j*mm_] -= (~rhs)(i,j);
            }
         }
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true>::subAssign( const SparseMatrix<MT,true>& rhs )
{
   for( size_t j=0UL; j<(~rhs).columns(); ++j )
      for( typename MT::ConstIterator element=(~rhs).begin(j); element!=(~rhs).end(j); ++element )
         v_[element->index()+j*mm_] -= element->value();
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
template< typename Type >  // Data type of the matrix
template< typename MT >    // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true>::subAssign( const SparseMatrix<MT,false>& rhs )
{
   for( size_t i=0UL; i<(~rhs).rows(); ++i )
      for( typename MT::ConstIterator element=(~rhs).begin(i); element!=(~rhs).end(i); ++element )
         v_[i+element->index()*mm_] -= element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  DYNAMICMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DynamicMatrix operators */
//@{
template< typename Type, bool SO >
inline void reset( DynamicMatrix<Type,SO>& m );

template< typename Type, bool SO >
inline void clear( DynamicMatrix<Type,SO>& m );

template< typename Type, bool SO >
inline bool isDefault( const DynamicMatrix<Type,SO>& m );

template< typename Type, bool SO >
inline void swap( DynamicMatrix<Type,SO>& a, DynamicMatrix<Type,SO>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dense matrix.
// \ingroup dynamic_matrix
//
// \param m The dense matrix to be resetted.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void reset( DynamicMatrix<Type,SO>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dense matrix.
// \ingroup dynamic_matrix
//
// \param m The dense matrix to be cleared.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void clear( DynamicMatrix<Type,SO>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense matrix is in default state.
// \ingroup dynamic_matrix
//
// \param m The dense matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
//
// This function checks whether the matrix is in default state. For instance, in case the
// matrix is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all matrix elements are 0 and \a false in case any matrix element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int> A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline bool isDefault( const DynamicMatrix<Type,SO>& m )
{
   if( SO == rowMajor ) {
      for( size_t i=0UL; i<m.rows(); ++i )
         for( size_t j=0UL; j<m.columns(); ++j )
            if( !isDefault( m(i,j) ) ) return false;
   }
   else {
      for( size_t j=0UL; j<m.columns(); ++j )
         for( size_t i=0UL; i<m.rows(); ++i )
            if( !isDefault( m(i,j) ) ) return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup dynamic_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void swap( DynamicMatrix<Type,SO>& a, DynamicMatrix<Type,SO>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO >
struct IsResizable< DynamicMatrix<T,SO> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};

template< typename T, bool SO >
struct IsResizable< const DynamicMatrix<T,SO> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};

template< typename T, bool SO >
struct IsResizable< volatile DynamicMatrix<T,SO> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};

template< typename T, bool SO >
struct IsResizable< const volatile DynamicMatrix<T,SO> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename T2, size_t M, size_t N >
struct AddTrait< DynamicMatrix<T1,SO>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticMatrix< typename AddTrait<T1,T2>::Type, M, N, SO >  Type;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct AddTrait< DynamicMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   typedef StaticMatrix< typename AddTrait<T1,T2>::Type, M, N, false >  Type;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct AddTrait< StaticMatrix<T1,M,N,SO>, DynamicMatrix<T2,SO> >
{
   typedef StaticMatrix< typename AddTrait<T1,T2>::Type, M, N, SO >  Type;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct AddTrait< StaticMatrix<T1,M,N,SO1>, DynamicMatrix<T2,SO2> >
{
   typedef StaticMatrix< typename AddTrait<T1,T2>::Type, M, N, false >  Type;
};

template< typename T1, bool SO, typename T2 >
struct AddTrait< DynamicMatrix<T1,SO>, DynamicMatrix<T2,SO> >
{
   typedef DynamicMatrix< typename AddTrait<T1,T2>::Type , SO >  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct AddTrait< DynamicMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   typedef DynamicMatrix< typename AddTrait<T1,T2>::Type , false >  Type;
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
template< typename T1, bool SO, typename T2, size_t M, size_t N >
struct SubTrait< DynamicMatrix<T1,SO>, StaticMatrix<T2,M,N,SO> >
{
   typedef StaticMatrix< typename SubTrait<T1,T2>::Type, M, N, SO >  Type;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct SubTrait< DynamicMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   typedef StaticMatrix< typename SubTrait<T1,T2>::Type, M, N, false >  Type;
};

template< typename T1, size_t M, size_t N, bool SO, typename T2 >
struct SubTrait< StaticMatrix<T1,M,N,SO>, DynamicMatrix<T2,SO> >
{
   typedef StaticMatrix< typename SubTrait<T1,T2>::Type, M, N, SO >  Type;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct SubTrait< StaticMatrix<T1,M,N,SO1>, DynamicMatrix<T2,SO2> >
{
   typedef StaticMatrix< typename SubTrait<T1,T2>::Type, M, N, false >  Type;
};

template< typename T1, bool SO, typename T2 >
struct SubTrait< DynamicMatrix<T1,SO>, DynamicMatrix<T2,SO> >
{
   typedef DynamicMatrix< typename SubTrait<T1,T2>::Type , SO >  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct SubTrait< DynamicMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   typedef DynamicMatrix< typename SubTrait<T1,T2>::Type , false >  Type;
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
template< typename T1, bool SO, typename T2 >
struct MultTrait< DynamicMatrix<T1,SO>, T2 >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, SO >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};

template< typename T1, typename T2, bool SO >
struct MultTrait< T1, DynamicMatrix<T2,SO> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, SO >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T1 );
};

template< typename T1, bool SO, typename T2, size_t N >
struct MultTrait< DynamicMatrix<T1,SO>, StaticVector<T2,N,false> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, size_t N, typename T2, bool SO >
struct MultTrait< StaticVector<T1,N,true>, DynamicMatrix<T2,SO> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, true >  Type;
};

template< typename T1, bool SO, typename T2, size_t N >
struct MultTrait< DynamicMatrix<T1,SO>, HybridVector<T2,N,false> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, size_t N, typename T2, bool SO >
struct MultTrait< HybridVector<T1,N,true>, DynamicMatrix<T2,SO> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, true >  Type;
};

template< typename T1, bool SO, typename T2 >
struct MultTrait< DynamicMatrix<T1,SO>, DynamicVector<T2,false> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< DynamicVector<T1,true>, DynamicMatrix<T2,SO> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, true >  Type;
};

template< typename T1, bool SO, typename T2 >
struct MultTrait< DynamicMatrix<T1,SO>, CompressedVector<T2,false> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< CompressedVector<T1,true>, DynamicMatrix<T2,SO> >
{
   typedef DynamicVector< typename MultTrait<T1,T2>::Type, true >  Type;
};

template< typename T1, bool SO1, typename T2, size_t M, size_t N, bool SO2 >
struct MultTrait< DynamicMatrix<T1,SO1>, StaticMatrix<T2,M,N,SO2> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, SO1 >  Type;
};

template< typename T1, size_t M, size_t N, bool SO1, typename T2, bool SO2 >
struct MultTrait< StaticMatrix<T1,M,N,SO1>, DynamicMatrix<T2,SO2> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, SO1 >  Type;
};

template< typename T1, bool SO1, typename T2, bool SO2 >
struct MultTrait< DynamicMatrix<T1,SO1>, DynamicMatrix<T2,SO2> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, SO1 >  Type;
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
template< typename T1, bool SO, typename T2 >
struct DivTrait< DynamicMatrix<T1,SO>, T2 >
{
   typedef DynamicMatrix< typename DivTrait<T1,T2>::Type , SO >  Type;
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
template< typename T1, bool SO, typename T2 >
struct MathTrait< DynamicMatrix<T1,SO>, DynamicMatrix<T2,SO> >
{
   typedef DynamicMatrix< typename MathTrait<T1,T2>::HighType, SO >  HighType;
   typedef DynamicMatrix< typename MathTrait<T1,T2>::LowType , SO >  LowType;
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
template< typename T1, bool SO >
struct SubmatrixTrait< DynamicMatrix<T1,SO> >
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
template< typename T1, bool SO >
struct RowTrait< DynamicMatrix<T1,SO> >
{
   typedef DynamicVector<T1,true>  Type;
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
template< typename T1, bool SO >
struct ColumnTrait< DynamicMatrix<T1,SO> >
{
   typedef DynamicVector<T1,false>  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
