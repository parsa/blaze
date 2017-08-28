//=================================================================================================
/*!
//  \file blaze/math/views/band/BaseTemplate.h
//  \brief Header file for the implementation of the Band base template
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_BAND_BASETEMPLATE_H_
#define _BLAZE_MATH_VIEWS_BAND_BASETEMPLATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
template< typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BIs >  // Band indices
class BandImpl
{};
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup band Band
// \ingroup views
*/
/*!\brief Reference to a specific band of a dense or sparse matrix.
// \ingroup band
//
// The Band template represents a reference to a specific band of a dense or sparse matrix
// primitive. The type of the matrix is specified via the first template parameter. Additionally
// it is possible to specify the band index at compile time as the second template parameter:

   \code
   template< typename MT, ptrdiff_t... BIs >
   class Band;
   \endcode

//  - MT: specifies the type of the matrix primitive. Band can be used with every matrix primitive,
//        but does not work with any matrix expression type.
//  - BIs: specifies zero or one band indices. This optional parameter can be used to provide the
//         band index at compile time.
//
//
// \n \section band_setup Setup of Bands
//
// \image html band.png
// \image latex band.eps "Band view" width=250pt
//
// A reference to a dense or sparse band can be created very conveniently via the \c band()
// function. The band index must be in the range from \f$[1-M..N-1]\f$, where \c M is the total
// number of rows and \c N is the total number of columns, and can be specified both at compile
// time or at runtime:

   \code
   using DenseMatrixType = blaze::DynamicMatrix<double,blaze::rowMajor>;

   DenseMatrixType A;
   // ... Resizing and initialization

   // Creating a reference to the 1st lower band of matrix A (compile time index)
   blaze::Band<DenseMatrixType,-1L> band1 = band<-1L>( A );

   // Creating a reference to the 2nd upper band of matrix A (runtime index)
   blaze::Band<DenseMatrixType> band2 = band( A, 2L );
   \endcode

// The resulting reference can be treated as any other vector, i.e. it can be assigned to, it can
// be copied from, and it can be used in arithmetic operations. By default, bands are considered
// column vectors, but this setting can be changed via the \c defaultTransposeFlag switch. The
// reference can also be used on both sides of an assignment: The band can either be used as an
// alias to grant write access to a specific band of a matrix primitive on the left-hand side of
// an assignment or to grant read-access to a specific band of a matrix primitive or expression
// on the right-hand side of an assignment. The following example demonstrates this in detail:

   \code
   using DenseVectorType  = blaze::DynamicVector<double,blaze::rowVector>;
   using SparseVectorType = blaze::CompressedVector<double,blaze::rowVector>;
   using DenseMatrixType  = blaze::DynamicMatrix<double,blaze::rowMajor>;
   using SparseMatrixType = blaze::CompressedMatrix<double,blaze::rowMajor>;

   DenseVectorType  x;
   SparseVectorType y;
   DenseMatrixType  A, B;
   SparseMatrixType C, D;
   // ... Resizing and initialization

   // Setting the 2nd upper band of matrix A to x
   blaze::Band<DenseMatrixType> band2 = band( A, 2L );
   band2 = x;

   // Setting the 3rd upper band of matrix B to y
   band( B, 3L ) = y;

   // Setting x to the 2nd lower band of the result of the matrix multiplication
   x = band( A * B, -2L );

   // Setting y to the 2nd upper band of the result of the sparse matrix multiplication
   y = band( C * D, 2L );
   \endcode

// The \c band() function can be used on any dense or sparse matrix, including expressions, as
// illustrated by the source code example. However, bands cannot be instantiated for expression
// types, but only for matrix primitives, respectively, i.e. for matrix types that offer write
// access.
//
//
// \n \section band_element_access Element access
//
// A dense or sparse band can be used like any other vector. For instance, the elements of a band
// can be directly accessed with the subscript operator:

   \code
   using MatrixType = blaze::DynamicMatrix<double,blaze::rowMajor>;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the 4th upper band of matrix A
   blaze::Band<MatrixType> band4 = band( A, 4L );

   // Setting the 1st element of the dense band, which corresponds
   // to the 1st element in the 4th upper band of matrix A
   band4[1] = 2.0;
   \endcode

// The numbering of the band elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of elements of the referenced band. Alternatively, the elements of a band
// can be traversed via iterators. Just as with vectors, in case of non-const band, \c begin() and
// \c end() return an Iterator, which allows a manipulation of the non-zero values, in case of
// constant band a ConstIterator is returned:

   \code
   using MatrixType = blaze::DynamicMatrix<int,blaze::rowMajor>;
   using BandType   = blaze::Band<MatrixType>;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 5th upper band of matrix A
   BandType band5 = band( A, 5L );

   for( BandType::Iterator it=band5.begin(); it!=band5.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense band value
      ... = *it;  // OK: Read access to the dense band value.
   }

   for( BandType::ConstIterator it=band5.begin(); it!=band5.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense band value.
   }
   \endcode

   \code
   using MatrixType = blaze::CompressedMatrix<int,blaze::rowMajor>;
   using BandType   = blaze::Band<MatrixType>;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 5th band of matrix A
   BandType band5 = band( A, 5L );

   for( BandType::Iterator it=band5.begin(); it!=band5.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( BandType::ConstIterator it=band5.begin(); it!=band5.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_band_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse band can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   using MatrixType = blaze::CompressedMatrix<double,blaze::rowMajor>;
   MatrixType A( 10UL, 100UL );  // Non-initialized 10x100 matrix

   using BandType = blaze::Band<MatrixType>;
   BandType diag( band( A, 0L ) );  // Reference to the diagonal of A

   // The subscript operator provides access to all possible elements of the sparse band,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse band, the element is inserted into the band.
   diag[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the band it is inserted into the band, if it is already contained in
   // the band its value is modified.
   diag.set( 45UL, -1.2 );

   // An alternative for inserting elements into the band is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the band.
   diag.insert( 50UL, 3.7 );
   \endcode

// \n \section band_common_operations Common Operations
//
// The current number of band elements can be obtained via the \c size() function, the current
// capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since bands are references to specific bands of a matrix,
// several operations are not possible, such as resizing and swapping. The following example
// shows this by means of a dense band view:

   \code
   using MatrixType = blaze::DynamicMatrix<int,blaze::rowMajor>;
   using BandType = blaze::Band<MatrixType>;

   MatrixType A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd upper band of matrix A
   BandType band2 = band( A, 2L );

   band2.size();          // Returns the number of elements in the band
   band2.capacity();      // Returns the capacity of the band
   band2.nonZeros();      // Returns the number of non-zero elements contained in the band

   band2.resize( 84UL );  // Compilation error: Cannot resize a single band of a matrix

   BandType band3 = band( A, 3L );
   swap( band2, band3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section band_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse bands can be used in all arithmetic operations that any other dense or
// sparse vector can be used in. The following example gives an impression of the use of dense
// bands within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse bands with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::columnVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::columnVector> c( 2UL );
   c[1] = 3.0;

   using DenseMatrix = blaze::DynamicMatrix<double,blaze::rowMajor>;
   DenseMatrix A( 4UL, 2UL );  // Non-initialized 4x2 matrix

   using BandType = blaze::Band<DenseMatrix>;
   BandType band1( band( A, 1L ) );  // Reference to the 1st upper band of A
   BandType diag ( band( A, 0L ) );  // Reference to the diagonal of A

   band1[0] = 0.0;      // Manual initialization of the 1st upper band of A
   diag = 1.0;          // Homogeneous initialization of the diagonal of A
   band( A, -1L ) = a;  // Dense vector initialization of the 1st lower band of A
   band( A, -2L ) = c;  // Sparse vector initialization of the 2nd lower band of A

   b = diag + a;               // Dense vector/dense vector addition
   b = c + band( A, -1L );     // Sparse vector/dense vector addition
   b = diag * band( A, -2L );  // Component-wise vector multiplication

   band( A, -1L ) *= 2.0;     // In-place scaling of the 1st upper band
   b = band( A, -1L ) * 2.0;  // Scaling of the 1st upper band
   b = 2.0 * band( A, -1L );  // Scaling of the 1st upper band

   band( A, -2L ) += a;              // Addition assignment
   band( A, -2L ) -= c;              // Subtraction assignment
   band( A, -2L ) *= band( A, 0L );  // Multiplication assignment

   double scalar = trans( c ) * band( A, -1L );  // Scalar/dot/inner product between two vectors

   A = band( A, -1L ) * trans( c );  // Outer product between two vectors
   \endcode
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
using Band = BandImpl< MT
                     , defaultTransposeFlag
                     , IsDenseMatrix<MT>::value
                     , IsMatMatMultExpr<MT>::value
                     , BIs... >;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reference to a specific band of a dense matrix.
// \ingroup band
//
// The DenseBand template represents a reference to a specific band of a dense matrix primitive.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
using DenseBand = BandImpl< MT
                          , defaultTransposeFlag
                          , true
                          , IsMatMatMultExpr<MT>::value
                          , BIs... >;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reference to a specific band of a sparse matrix.
// \ingroup band
//
// The SparseBand template represents a reference to a specific band of a sparse matrix primitive.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
using SparseBand = BandImpl< MT
                           , defaultTransposeFlag
                           , false
                           , IsMatMatMultExpr<MT>::value
                           , BIs... >;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\defgroup diagonal Diagonal
// \ingroup views
*/
/*!\brief Reference to the diagonal of a dense or sparse matrix.
// \ingroup diagonal
//
// The Diagonal template represents a reference to a specific diagonal of a dense or sparse matrix
// primitive. The type of the matrix is specified via template parameter:

   \code
   template< typename MT >
   class Diagonal;
   \endcode

//  - MT: specifies the type of the matrix primitive. Diagonal can be used with every matrix
//        primitive, but does not work with any matrix expression type.
//
//
// \n \section diagonal_setup Setup of Diagonals
//
// A reference to a dense or sparse diagonal can be created very conveniently via the \c diagonal()
// function:

   \code
   using DenseMatrixType = blaze::DynamicMatrix<double,blaze::rowMajor>;

   DenseMatrixType A;
   // ... Resizing and initialization

   // Creating a reference to the diagonal of matrix A
   blaze::Diagonal<DenseMatrixType> diag = diagonal( A );
   \endcode

// The resulting reference can be treated as any other vector, i.e. it can be assigned to, it
// can be copied from, and it can be used in arithmetic operations. By default, diagonals are
// considered column vectors, but this setting can be changed via the \c defaultTransposeFlag
// switch. The reference can also be used on both sides of an assignment: The diagonal can either
// be used as an alias to grant write access to a specific diagonal of a matrix primitive on the
// left-hand side of an assignment or to grant read-access to a specific diagonal of a matrix
// primitive or expression on the right-hand side of an assignment. The following example
// demonstrates this in detail:

   \code
   using DenseVectorType  = blaze::DynamicVector<double,blaze::rowVector>;
   using SparseVectorType = blaze::CompressedVector<double,blaze::rowVector>;
   using DenseMatrixType  = blaze::DynamicMatrix<double,blaze::rowMajor>;
   using SparseMatrixType = blaze::CompressedMatrix<double,blaze::rowMajor>;

   DenseVectorType  x;
   SparseVectorType y;
   DenseMatrixType  A, B;
   SparseMatrixType C, D;
   // ... Resizing and initialization

   // Setting the diagonal of matrix A to x
   blaze::Diagonal<DenseMatrixType> diag = diagonal( A );
   diag = x;

   // Setting the diagonal of matrix B to y
   diagonal( B ) = y;

   // Setting x to the diagonal of the result of the matrix multiplication
   x = diagonal( A * B );

   // Setting y to the digaonal of the result of the sparse matrix multiplication
   y = diagonal( C * D );
   \endcode

// The \c diagaonl() function can be used on any dense or sparse matrix, including expressions, as
// illustrated by the source code example. However, diagonals cannot be instantiated for expression
// types, but only for matrix primitives, respectively, i.e. for matrix types that offer write
// access.
//
//
// \n \section diagonal_element_access Element access
//
// A dense or sparse diagonal can be used like any other vector. For instance, the elements of a
// diagonal can be directly accessed with the subscript operator:

   \code
   using MatrixType = blaze::DynamicMatrix<double,blaze::rowMajor>;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the diagonal of matrix A
   blaze::Diagonal<MatrixType> diag = diagonal( A );

   // Setting the 1st element of the dense diagonal, which corresponds
   // to the 1st element on the diagonal of matrix A
   diag[1] = 2.0;
   \endcode

// The numbering of the diagonal elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of elements of the referenced diagonal. Alternatively, the elements of a
// diagonal can be traversed via iterators. Just as with vectors, in case of non-const diagonal,
// \c begin() and \c end() return an Iterator, which allows a manipulation of the non-zero values,
// in case of constant diagonal a ConstIterator is returned:

   \code
   using MatrixType   = blaze::DynamicMatrix<int,blaze::rowMajor>;
   using DiagonalType = blaze::Diagonal<MatrixType>;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the diagonal of matrix A
   DiagonalType diag = Diagonal( A );

   for( DiagonalType::Iterator it=diag.begin(); it!=diag.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense diagonal value
      ... = *it;  // OK: Read access to the dense diagonal value.
   }

   for( DiagonalType::ConstIterator it=diag.begin(); it!=diag.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense diagonal value.
   }
   \endcode

   \code
   using MatrixType   = blaze::CompressedMatrix<int,blaze::rowMajor>;
   using DiagonalType = blaze::Diagonal<MatrixType>;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the diagonal of matrix A
   DiagonalType diag = diagonal( A );

   for( DiagonalType::Iterator it=diag.begin(); it!=diag.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( DiagonalType::ConstIterator it=diag.begin(); it!=diag.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section sparse_diagonal_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse diagonal can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   using MatrixType = blaze::CompressedMatrix<double,blaze::rowMajor>;
   MatrixType A( 10UL, 100UL );  // Non-initialized 10x100 matrix

   using DiagonalType = blaze::Diaganal<MatrixType>;
   DiagonalType diag( diagonal( A ) );  // Reference to the diagonal of A

   // The subscript operator provides access to all possible elements of the sparse diagonal,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse diagonal, the element is inserted into the
   // diagonal.
   diag[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element is
   // not contained in the diagonal it is inserted into the diagonal, if it is already contained
   // in the diagonal its value is modified.
   diag.set( 45UL, -1.2 );

   // An alternative for inserting elements into the diagonal is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the diagonal.
   diag.insert( 50UL, 3.7 );
   \endcode

// \n \section diagonal_common_operations Common Operations
//
// The current number of diagonal elements can be obtained via the \c size() function, the
// current capacity via the \c capacity() function, and the number of non-zero elements via the
// \c nonZeros() function. However, since diagonals are references to specific diagonals of a
// matrix, several operations are not possible, such as resizing and swapping. The following
// example shows this by means of a dense diagonal view:

   \code
   using MatrixType   = blaze::DynamicMatrix<int,blaze::rowMajor>;
   using DiagonalType = blaze::Diagonal<MatrixType>;

   MatrixType A( 42UL, 42UL );
   MatrixType B( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the diagonal of matrix A
   DiagonalType diag1 = diagonal( A );

   diag1.size();          // Returns the number of elements in the diagonal
   diag1.capacity();      // Returns the capacity of the diagonal
   diag1.nonZeros();      // Returns the number of non-zero elements contained in the diagonal

   diag1.resize( 84UL );  // Compilation error: Cannot resize the diagonal of a matrix

   DiagonalType diag2 = diagonal( B );
   swap( diag1, diag2 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section diagonal_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse diagonals can be used in all arithmetic operations that any other dense
// or sparse vector can be used in. The following example gives an impression of the use of dense
// diagonals within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse diagonals with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::columnVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::columnVector> c( 2UL );
   c[1] = 3.0;

   using DenseMatrix = blaze::DynamicMatrix<double,blaze::rowMajor>;
   DenseMatrix A( 4UL, 2UL );  // Non-initialized 4x2 matrix
   DenseMatrix B( 2UL, 4UL );  // Non-initialized 2x4 matrix

   using DiagonalType = blaze::Diagonal<DenseMatrix>;
   DiagonalType diag1( diagonal( A ) );  // Reference to the diagonal of A
   DiagonalType diag2( diagonal( B ) );  // Reference to the diagonal of B

   diag1[0] = 0.0;     // Manual initialization of the diagonal of A
   diag2 = 1.0;        // Homogeneous initialization of the diagonal of A
   diagonal( A ) = a;  // Dense vector initialization of the diagonal of A
   diagonal( A ) = c;  // Sparse vector initialization of the diagonal of A

   b = diag2 + a;              // Dense vector/dense vector addition
   b = c + diagonal( A );      // Sparse vector/dense vector addition
   b = diag2 * diagonal( A );  // Component-wise vector multiplication

   diagonal( A ) *= 2.0;     // In-place scaling of the diagonal
   b = diagonal( A ) * 2.0;  // Scaling of the diagonal
   b = 2.0 * diagonal( A );  // Scaling of the diagonal

   diagonal( A ) += a;              // Addition assignment
   diagonal( A ) -= c;              // Subtraction assignment
   diagonal( A ) *= diagonal( A );  // Multiplication assignment

   double scalar = trans( c ) * diagonal( A );  // Scalar/dot/inner product between two vectors

   A = diagonal( A ) * trans( c );  // Outer product between two vectors
   \endcode
*/
template< typename MT >
using Diagonal = Band<MT,0L>;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reference to a specific diagonal of a dense matrix.
// \ingroup diagonal
//
// The DenseDiagonal template represents a reference to a specific diagonal of a dense matrix
// primitive.
*/
template< typename MT >   // Type of the matrix
using DenseDiagonal = DenseBand<MT,0L>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reference to a specific diagonal of a sparse matrix.
// \ingroup diagonal
//
// The SparseDiagonal template represents a reference to a specific diagonal of a sparse matrix
// primitive.
*/
template< typename MT >  // Type of the matrix
using SparseDiagonal = SparseBand<MT,0L>;
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
