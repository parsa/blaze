//=================================================================================================
/*!
//  \file blaze/math/views/Submatrix.h
//  \brief Header file for the implementation of the Submatrix view
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

#ifndef _BLAZE_MATH_VIEWS_SUBMATRIX_H_
#define _BLAZE_MATH_VIEWS_SUBMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/expressions/DeclExpr.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnsTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/RowsTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/column/ColumnData.h>
#include <blaze/math/views/row/RowData.h>
#include <blaze/math/views/submatrix/BaseTemplate.h>
#include <blaze/math/views/submatrix/Dense.h>
#include <blaze/math/views/submatrix/Sparse.h>
#include <blaze/math/views/Subvector.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/TypeList.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup submatrix
//
// \param matrix The matrix containing the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given matrix.
// The following example demonstrates the creation of a dense and sparse submatrix:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D;
   blaze::CompressedMatrix<int,blaze::columnMajor> S;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix<0UL,16UL,8UL,4UL>( D );

   // Creating a sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix<2UL,4UL,7UL,3UL>( S );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix<0UL,16UL,8UL,4UL>( D, unchecked );
   auto ssm = submatrix<2UL,4UL,7UL,3UL>( S, unchecked );
   \endcode

// Please note that this function creates an unaligned dense or sparse submatrix. For instance,
// the creation of the dense submatrix is equivalent to the following function call:

   \code
   auto dsm = submatrix<unaligned,0UL,16UL,8UL,4UL>( D );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions. However, especially in case of dense submatrices this may
// result in considerable performance improvements. In order to create an aligned submatrix the
// following function call has to be used:

   \code
   auto dsm = submatrix<aligned,0UL,16UL,8UL,4UL>( D );
   \endcode

// Note however that in this case the given compile time arguments \a I, \a J, \a M, and \a N are
// subject to additional checks to guarantee proper alignment.
*/
template< size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( Matrix<MT,SO>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned,I,J,M,N>( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given constant matrix.
// \ingroup submatrix
//
// \param matrix The constant matrix containing the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given constant
// matrix. The following example demonstrates the creation of a dense and sparse submatrix:

   \code
   const blaze::DynamicMatrix<double,blaze::rowMajor> D( ... );
   const blaze::CompressedMatrix<int,blaze::columnMajor> S( ... );

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix<0UL,16UL,8UL,4UL>( D );

   // Creating a sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix<2UL,4UL,7UL,3UL>( S );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix<0UL,16UL,8UL,4UL>( D, unchecked );
   auto ssm = submatrix<2UL,4UL,7UL,3UL>( S, unchecked );
   \endcode

// Please note that this function creates an unaligned dense or sparse submatrix. For instance,
// the creation of the dense submatrix is equivalent to the following three function calls:

   \code
   auto dsm = submatrix<unaligned,0UL,16UL,8UL,4UL>( D );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions. However, especially in case of dense submatrices this may
// result in considerable performance improvements. In order to create an aligned submatrix the
// following function call has to be used:

   \code
   auto dsm = submatrix<aligned,0UL,16UL,8UL,4UL>( D );
   \endcode

// Note however that in this case the given compile time arguments \a I, \a J, \a M, and \a N are
// subject to additional checks to guarantee proper alignment.
*/
template< size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto) submatrix( const Matrix<MT,SO>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned,I,J,M,N>( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given temporary matrix.
// \ingroup submatrix
//
// \param matrix The temporary matrix containing the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given
// temporary matrix. In case the submatrix is not properly specified (i.e. if the specified
// row or column is greater than the total number of rows or columns of the given matrix or
// the submatrix is specified beyond the number of rows or columns of the matrix) a
// \a std::invalid_argument exception is thrown.
*/
template< size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto) submatrix( Matrix<MT,SO>&& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned,I,J,M,N>( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup submatrix
//
// \param matrix The matrix containing the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given dense or sparse matrix, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned submatrix:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D;
   blaze::CompressedMatrix<int,blaze::columnMajor> S;
   // ... Resizing and initialization

   // Creating an aligned dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix<aligned,0UL,16UL,8UL,4UL>( D );

   // Creating an unaligned sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix<unaligned,2UL,4UL,7UL,3UL>( S );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix<aligned,0UL,16UL,8UL,4UL>( D, unchecked );
   auto ssm = submatrix<unaligned,2UL,4UL,7UL,3UL>( S, unchecked );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices
// pose additional alignment restrictions and the given \a I, and \a J arguments are subject
// to additional checks to guarantee proper alignment. However, especially in case of dense
// submatrices this may result in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of each row/column of the submatrix must be aligned. The following source code
// gives some examples for a double precision row-major dynamic matrix, assuming that padding is
// enabled and that AVX is available, which packs 4 \c double values into a SIMD vector:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0), i.e. the first element of each row is aligned (due to padding)
   auto dsm1 = submatrix<aligned,0UL,0UL,7UL,11UL>( D );

   // OK: First column is a multiple of 4, i.e. the first element of each row is aligned (due to padding)
   auto dsm2 = submatrix<aligned,3UL,12UL,8UL,16UL>( D );

   // OK: First column is a multiple of 4 and the submatrix includes the last row and column
   auto dsm3 = submatrix<aligned,4UL,0UL,9UL,17UL>( D );

   // Error: First column is not a multiple of 4, i.e. the first element is not aligned
   auto dsm4 = submatrix<aligned,2UL,3UL,12UL,12UL>( D );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto) submatrix( Matrix<MT,SO>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Submatrix_<MT,AF,I,J,M,N>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given constant matrix.
// \ingroup submatrix
//
// \param matrix The constant matrix containing the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given constant dense or sparse matrix, based on the specified alignment flag \a AF. The
// following example demonstrates the creation of both an aligned and unaligned submatrix:

   \code
   const blaze::DynamicMatrix<double,blaze::rowMajor> D( ... );
   const blaze::CompressedMatrix<int,blaze::columnMajor> S( ... );

   // Creating an aligned dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix<aligned,0UL,16UL,8UL,4UL>( D );

   // Creating an unaligned sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix<unaligned,2UL,4UL,7UL,3UL>( S );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix<aligned,0UL,16UL,8UL,4UL>( D, unchecked );
   auto ssm = submatrix<unaligned,2UL,4UL,7UL,3UL>( S, unchecked );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices
// pose additional alignment restrictions and the given \a I, and \a J arguments are subject
// to additional checks to guarantee proper alignment. However, especially in case of dense
// submatrices this may result in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of each row/column of the submatrix must be aligned. The following source code
// gives some examples for a double precision row-major dynamic matrix, assuming that padding is
// enabled and that AVX is available, which packs 4 \c double values into a SIMD vector:

   \code
   const blaze::DynamicMatrix<double,blaze::rowMajor> D( ... );

   // OK: Starts at position (0,0), i.e. the first element of each row is aligned (due to padding)
   auto dsm1 = submatrix<aligned,0UL,0UL,7UL,11UL>( D );

   // OK: First column is a multiple of 4, i.e. the first element of each row is aligned (due to padding)
   auto dsm2 = submatrix<aligned,3UL,12UL,8UL,16UL>( D );

   // OK: First column is a multiple of 4 and the submatrix includes the last row and column
   auto dsm3 = submatrix<aligned,4UL,0UL,9UL,17UL>( D );

   // Error: First column is not a multiple of 4, i.e. the first element is not aligned
   auto dsm4 = submatrix<aligned,2UL,3UL,12UL,12UL>( D );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto) submatrix( const Matrix<MT,SO>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Submatrix_<const MT,AF,I,J,M,N>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given temporary matrix.
// \ingroup submatrix
//
// \param matrix The temporary matrix containing the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given temporary dense or sparse matrix, based on the specified alignment flag \a AF. In
// case the submatrix is not properly specified (i.e. if the specified row or column is larger
// than the total number of rows or columns of the given matrix or the submatrix is specified
// beyond the number of rows or columns of the matrix) or any alignment restrictions are
// violated, a \a std::invalid_argument exception is thrown.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto) submatrix( Matrix<MT,SO>&& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Submatrix_<MT,AF,I,J,M,N>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup submatrix
//
// \param matrix The matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given matrix.
// The following example demonstrates the creation of a dense and sparse submatrix:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D;
   blaze::CompressedMatrix<int,blaze::columnMajor> S;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix( D, 0UL, 16UL, 8UL, 4UL );

   // Creating a sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix( D, 0UL, 16UL, 8UL, 4UL, unchecked );
   auto ssm = submatrix( S, 2UL, 4UL, 7UL, 3UL, unchecked );
   \endcode

// Please note that this function creates an unaligned dense or sparse submatrix. For instance,
// the creation of the dense submatrix is equivalent to the following function call:

   \code
   unaligned dsm = submatrix<unaligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions. However, especially in case of dense submatrices this may
// result in considerable performance improvements. In order to create an aligned submatrix the
// following function call has to be used:

   \code
   auto dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// Note however that in this case the given arguments \a row, \a column, \a m, and \a n are
// subject to additional checks to guarantee proper alignment.
*/
template< typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto)
   submatrix( Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned>( ~matrix, row, column, m, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given constant matrix.
// \ingroup submatrix
//
// \param matrix The constant matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given constant
// matrix. The following example demonstrates the creation of a dense and sparse submatrix:

   \code
   const blaze::DynamicMatrix<double,blaze::rowMajor> D( ... );
   const blaze::CompressedMatrix<int,blaze::columnMajor> S( ... );

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix( D, 0UL, 16UL, 8UL, 4UL );

   // Creating a sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix( D, 0UL, 16UL, 8UL, 4UL, unchecked );
   auto ssm = submatrix( S, 2UL, 4UL, 7UL, 3UL, unchecked );
   \endcode

// Please note that this function creates an unaligned dense or sparse submatrix. For instance,
// the creation of the dense submatrix is equivalent to the following three function calls:

   \code
   auto dsm = submatrix<unaligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions. However, especially in case of dense submatrices this may
// result in considerable performance improvements. In order to create an aligned submatrix the
// following function call has to be used:

   \code
   auto dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// Note however that in this case the given arguments \a row, \a column, \a m, and \a n are
// subject to additional checks to guarantee proper alignment.
*/
template< typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto)
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned>( ~matrix, row, column, m, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given temporary matrix.
// \ingroup submatrix
//
// \param matrix The temporary matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given
// temporary matrix. In case the submatrix is not properly specified (i.e. if the specified
// row or column is greater than the total number of rows or columns of the given matrix or
// the submatrix is specified beyond the number of rows or columns of the matrix) a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto)
   submatrix( Matrix<MT,SO>&& matrix, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned>( ~matrix, row, column, m, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup submatrix
//
// \param matrix The matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given dense or sparse matrix, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned submatrix:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D;
   blaze::CompressedMatrix<int,blaze::columnMajor> S;
   // ... Resizing and initialization

   // Creating an aligned dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );

   // Creating an unaligned sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix<unaligned>( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL, unchecked );
   auto ssm = submatrix<unaligned>( S, 2UL, 4UL, 7UL, 3UL, unchecked );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions and the given \a row, and \a column arguments are subject to
// additional checks to guarantee proper alignment. However, especially in case of dense submatrices
// this may result in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of each row/column of the submatrix must be aligned. The following source code
// gives some examples for a double precision row-major dynamic matrix, assuming that padding is
// enabled and that AVX is available, which packs 4 \c double values into a SIMD vector:

   \code
   blaze::DynamicMatrix<double,blaze::rowMajor> D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0), i.e. the first element of each row is aligned (due to padding)
   auto dsm1 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 11UL );

   // OK: First column is a multiple of 4, i.e. the first element of each row is aligned (due to padding)
   auto dsm2 = submatrix<aligned>( D, 3UL, 12UL, 8UL, 16UL );

   // OK: First column is a multiple of 4 and the submatrix includes the last row and column
   auto dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First column is not a multiple of 4, i.e. the first element is not aligned
   auto dsm4 = submatrix<aligned>( D, 2UL, 3UL, 12UL, 12UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto)
   submatrix( Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Submatrix_<MT,AF>;
   return ReturnType( ~matrix, row, column, m, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given constant matrix.
// \ingroup submatrix
//
// \param matrix The constant matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given dense or sparse matrix, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned submatrix:

   \code
   const blaze::DynamicMatrix<double,blaze::rowMajor> D( ... );
   const blaze::CompressedMatrix<int,blaze::columnMajor> S( ... );

   // Creating an aligned dense submatrix of size 8x4, starting in row 0 and column 16
   auto dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );

   // Creating an unaligned sparse submatrix of size 7x3, starting in row 2 and column 4
   auto ssm = submatrix<unaligned>( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// By default, the provided submatrix arguments are checked at runtime. In case the submatrix
// is not properly specified (i.e. if the specified row or column is larger than the total number
// of rows or columns of the given matrix or the submatrix is specified beyond the number of rows
// or columns of the matrix) a \a std::invalid_argument exception is thrown. The checks can be
// skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL, unchecked );
   auto ssm = submatrix<unaligned>( S, 2UL, 4UL, 7UL, 3UL, unchecked );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions and the given \a row, and \a column arguments are subject to
// additional checks to guarantee proper alignment. However, especially in case of dense submatrices
// this may result in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of each row/column of the submatrix must be aligned. The following source code
// gives some examples for a double precision row-major dynamic matrix, assuming that padding is
// enabled and that AVX is available, which packs 4 \c double values into a SIMD vector:

   \code
   const blaze::DynamicMatrix<double,blaze::rowMajor> D( ... );

   // OK: Starts at position (0,0), i.e. the first element of each row is aligned (due to padding)
   auto dsm1 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 11UL );

   // OK: First column is a multiple of 4, i.e. the first element of each row is aligned (due to padding)
   auto dsm2 = submatrix<aligned>( D, 3UL, 12UL, 8UL, 16UL );

   // OK: First column is a multiple of 4 and the submatrix includes the last row and column
   auto dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First column is not a multiple of 4, i.e. the first element is not aligned
   auto dsm4 = submatrix<aligned>( D, 2UL, 3UL, 12UL, 12UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto)
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Submatrix_<const MT,AF>;
   return ReturnType( ~matrix, row, column, m, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given temporary matrix.
// \ingroup submatrix
//
// \param matrix The temporary matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args Optional submatrix arguments.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given temporary dense or sparse matrix, based on the specified alignment flag \a AF. In
// case the submatrix is not properly specified (i.e. if the specified row or column is larger
// than the total number of rows or columns of the given matrix or the submatrix is specified
// beyond the number of rows or columns of the matrix) or any alignment restrictions are
// violated, a \a std::invalid_argument exception is thrown.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , typename... RSAs >  // Option submatrix arguments
inline decltype(auto)
   submatrix( Matrix<MT,SO>&& matrix, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Submatrix_<MT,AF>;
   return ReturnType( ~matrix, row, column, m, n, args... );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix addition.
// \ingroup submatrix
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime submatrix arguments
// \return View on the specified submatrix of the addition.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix addition.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatMatAddExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF,CSAs...>( (~matrix).leftOperand(), args... ) +
          submatrix<AF,CSAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix subtraction.
// \ingroup submatrix
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime submatrix arguments
// \return View on the specified submatrix of the subtraction.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix subtraction.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatMatSubExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF,CSAs...>( (~matrix).leftOperand(), args... ) -
          submatrix<AF,CSAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given Schur product.
// \ingroup submatrix
//
// \param matrix The constant Schur product.
// \param args The runtime submatrix arguments
// \return View on the specified submatrix of the Schur product.
//
// This function returns an expression representing the specified submatrix of the given Schur
// product.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const SchurExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF,CSAs...>( (~matrix).leftOperand(), args... ) %
          submatrix<AF,CSAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix multiplication.
// \ingroup submatrix
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the multiplication.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix multiplication.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatMatMultExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using MT1 = RemoveReference_< LeftOperand_< MatrixType_<MT> > >;
   using MT2 = RemoveReference_< RightOperand_< MatrixType_<MT> > >;

   const SubmatrixData<CSAs...> sd( args... );

   decltype(auto) left ( (~matrix).leftOperand()  );
   decltype(auto) right( (~matrix).rightOperand() );

   const size_t begin( max( ( IsUpper<MT1>::value )
                            ?( ( !AF && IsStrictlyUpper<MT1>::value )
                               ?( sd.row() + 1UL )
                               :( sd.row() ) )
                            :( 0UL )
                          , ( IsLower<MT2>::value )
                            ?( ( !AF && IsStrictlyLower<MT2>::value )
                               ?( sd.column() + 1UL )
                               :( sd.column() ) )
                            :( 0UL ) ) );
   const size_t end( min( ( IsLower<MT1>::value )
                          ?( ( IsStrictlyLower<MT1>::value && sd.rows() > 0UL )
                             ?( sd.row() + sd.rows() - 1UL )
                             :( sd.row() + sd.rows() ) )
                          :( left.columns() )
                        , ( IsUpper<MT2>::value )
                          ?( ( IsStrictlyUpper<MT2>::value && sd.columns() > 0UL )
                             ?( sd.column() + sd.columns() - 1UL )
                             :( sd.column() + sd.columns() ) )
                          :( left.columns() ) ) );

   const size_t diff( ( begin < end )?( end - begin ):( 0UL ) );

   return submatrix<AF>( left, sd.row(), begin, sd.rows(), diff ) *
          submatrix<AF>( right, begin, sd.column(), diff, sd.columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given outer product.
// \ingroup submatrix
//
// \param matrix The constant outer product.
// \return View on the specified submatrix of the outer product.
//
// This function returns an expression representing the specified submatrix of the given
// outer product.
*/
template< AlignmentFlag AF  // Alignment flag
        , size_t I          // Index of the first row
        , size_t J          // Index of the first column
        , size_t M          // Number of rows
        , size_t N          // Number of columns
        , typename MT >     // Matrix base type of the expression
inline decltype(auto) submatrix( const VecTVecMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,M>( (~matrix).leftOperand() ) *
          subvector<AF,J,N>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given outer product.
// \ingroup submatrix
//
// \param matrix The constant outer product.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the outer product.
//
// This function returns an expression representing the specified submatrix of the given
// outer product.
*/
template< AlignmentFlag AF  // Alignment flag
        , typename MT >     // Matrix base type of the expression
inline decltype(auto)
   submatrix( const VecTVecMultExpr<MT>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~matrix).leftOperand(), row, m ) *
          subvector<AF>( (~matrix).rightOperand(), column, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/scalar multiplication.
// \ingroup submatrix
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the multiplication.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/scalar multiplication.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatScalarMultExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF,CSAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/scalar division.
// \ingroup submatrix
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the division.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/scalar division.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatScalarDivExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF,CSAs...>( (~matrix).leftOperand(), args... ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given unary matrix map operation.
// \ingroup submatrix
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the unary map operation.
//
// This function returns an expression representing the specified submatrix of the given unary
// matrix map operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatMapExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( submatrix<AF,CSAs...>( (~matrix).operand(), args... ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given binary matrix map operation.
// \ingroup submatrix
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the binary map operation.
//
// This function returns an expression representing the specified submatrix of the given binary
// matrix map operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatMatMapExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( submatrix<AF,CSAs...>( (~matrix).leftOperand(), args... ),
               submatrix<AF,CSAs...>( (~matrix).rightOperand(), args... ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix evaluation operation.
// \ingroup submatrix
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the evaluation operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// evaluation operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatEvalExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( submatrix<AF,CSAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix serialization operation.
// \ingroup submatrix
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the serialization operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// serialization operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const MatSerialExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( submatrix<AF,CSAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix declaration operation.
// \ingroup submatrix
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime submatrix arguments.
// \return View on the specified submatrix of the declaration operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// declaration operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Runtime submatrix arguments
inline decltype(auto) submatrix( const DeclExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF,CSAs...>( (~matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix transpose operation.
// \ingroup submatrix
//
// \param matrix The constant matrix transpose operation.
// \param args Optional submatrix arguments.
// \return View on the specified submatrix of the transpose operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// transpose operation.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename MT         // Matrix base type of the expression
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( const MatTransExpr<MT>& matrix, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return trans( submatrix<AF,CSAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another submatrix.
// \ingroup submatrix
//
// \param sm The given submatrix
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , size_t I1           // Required index of the first row
        , size_t J1           // Required index of the first column
        , size_t M1           // Required number of rows
        , size_t N1           // Required number of columns
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J2           // Present index of the first column
        , size_t M2           // Present number of rows
        , size_t N2           // Present number of columns
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( Submatrix<MT,AF2,SO,DF,I2,J2,M2,N2>& sm, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_STATIC_ASSERT_MSG( I1 + M1 <= M2, "Invalid submatrix specification" );
   BLAZE_STATIC_ASSERT_MSG( J1 + N1 <= N2, "Invalid submatrix specification" );

   return submatrix<AF1,I1+I2,J1+J2,M1,N1>( sm.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another constant submatrix.
// \ingroup submatrix
//
// \param sm The given constant submatrix
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given constant
// submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , size_t I1           // Required index of the first row
        , size_t J1           // Required index of the first column
        , size_t M1           // Required number of rows
        , size_t N1           // Required number of columns
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J2           // Present index of the first column
        , size_t M2           // Present number of rows
        , size_t N2           // Present number of columns
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( const Submatrix<MT,AF2,SO,DF,I2,J2,M2,N2>& sm, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_STATIC_ASSERT_MSG( I1 + M1 <= M2, "Invalid submatrix specification" );
   BLAZE_STATIC_ASSERT_MSG( J1 + N1 <= N2, "Invalid submatrix specification" );

   return submatrix<AF1,I1+I2,J1+J2,M1,N1>( sm.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another temporary submatrix.
// \ingroup submatrix
//
// \param sm The given temporary submatrix
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given temporary
// submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , size_t I1           // Required index of the first row
        , size_t J1           // Required index of the first column
        , size_t M1           // Required number of rows
        , size_t N1           // Required number of columns
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J2           // Present index of the first column
        , size_t M2           // Present number of rows
        , size_t N2           // Present number of columns
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( Submatrix<MT,AF2,SO,DF,I2,J2,M2,N2>&& sm, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_STATIC_ASSERT_MSG( I1 + M1 <= M2, "Invalid submatrix specification" );
   BLAZE_STATIC_ASSERT_MSG( J1 + N1 <= N2, "Invalid submatrix specification" );

   return submatrix<AF1,I1+I2,J1+J2,M1,N1>( sm.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another submatrix.
// \ingroup submatrix
//
// \param sm The given submatrix
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( Submatrix<MT,AF2,SO,DF>& sm, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains< TypeList<RSAs...>, Unchecked >::value );

   if( isChecked ) {
      if( ( I + M > sm.rows() ) || ( J + N > sm.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I + M <= sm.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( J + N <= sm.columns(), "Invalid submatrix specification" );
   }

   return submatrix<AF1>( sm.operand(), sm.row() + I, sm.column() + J, M, N, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another constant submatrix.
// \ingroup submatrix
//
// \param sm The constant submatrix
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given constant
// submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( const Submatrix<MT,AF2,SO,DF>& sm, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains< TypeList<RSAs...>, Unchecked >::value );

   if( isChecked ) {
      if( ( I + M > sm.rows() ) || ( J + N > sm.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I + M <= sm.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( J + N <= sm.columns(), "Invalid submatrix specification" );
   }

   return submatrix<AF1>( sm.operand(), sm.row() + I, sm.column() + J, M, N, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another temporary submatrix.
// \ingroup submatrix
//
// \param sm The temporary submatrix
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given temporary
// submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto) submatrix( Submatrix<MT,AF2,SO,DF>&& sm, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains< TypeList<RSAs...>, Unchecked >::value );

   if( isChecked ) {
      if( ( I + M > sm.rows() ) || ( J + N > sm.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I + M <= sm.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( J + N <= sm.columns(), "Invalid submatrix specification" );
   }

   return submatrix<AF1>( sm.operand(), sm.row() + I, sm.column() + J, M, N, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another submatrix.
// \ingroup submatrix
//
// \param sm The given submatrix
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto)
   submatrix( Submatrix<MT,AF2,SO,DF,CSAs...>& sm, size_t row, size_t column,
              size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains< TypeList<RSAs...>, Unchecked >::value );

   if( isChecked ) {
      if( ( row + m > sm.rows() ) || ( column + n > sm.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row    + m <= sm.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column + n <= sm.columns(), "Invalid submatrix specification" );
   }

   return submatrix<AF1>( sm.operand(), sm.row() + row, sm.column() + column, m, n, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another constant submatrix.
// \ingroup submatrix
//
// \param sm The given constant submatrix
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given constant
// submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto)
   submatrix( const Submatrix<MT,AF2,SO,DF,CSAs...>& sm, size_t row, size_t column,
              size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains< TypeList<RSAs...>, Unchecked >::value );

   if( isChecked ) {
      if( ( row + m > sm.rows() ) || ( column + n > sm.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row    + m <= sm.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column + n <= sm.columns(), "Invalid submatrix specification" );
   }

   return submatrix<AF1>( sm.operand(), sm.row() + row, sm.column() + column, m, n, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of another temporary submatrix.
// \ingroup submatrix
//
// \param sm The given temporary submatrix
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the other submatrix.
//
// This function returns an expression representing the specified submatrix of the given temporary
// submatrix.
*/
template< AlignmentFlag AF1   // Required alignment flag
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF2   // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RSAs >  // Optional submatrix arguments
inline decltype(auto)
   submatrix( Submatrix<MT,AF2,SO,DF,CSAs...>&& sm, size_t row, size_t column,
              size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains< TypeList<RSAs...>, Unchecked >::value );

   if( isChecked ) {
      if( ( row + m > sm.rows() ) || ( column + n > sm.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row    + m <= sm.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column + n <= sm.columns(), "Invalid submatrix specification" );
   }

   return submatrix<AF1>( sm.operand(), sm.row() + row, sm.column() + column, m, n, args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (SUBVECTOR)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given matrix/vector multiplication.
// \ingroup submatrix
//
// \param vector The constant matrix/vector multiplication.
// \param args The runtime subvector arguments.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// matrix/vector multiplication.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time subvector arguments
        , typename VT         // Vector base type of the expression
        , typename... RSAs >  // Runtime subvector arguments
inline decltype(auto) subvector( const MatVecMultExpr<VT>& vector, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using MT = RemoveReference_< LeftOperand_< VectorType_<VT> > >;

   const SubvectorData<CSAs...> sd( args... );

   decltype(auto) left ( (~vector).leftOperand()  );
   decltype(auto) right( (~vector).rightOperand() );

   const size_t column( ( IsUpper<MT>::value )
                        ?( ( !AF && IsStrictlyUpper<MT>::value )?( sd.offset() + 1UL ):( sd.offset() ) )
                        :( 0UL ) );
   const size_t n( ( IsLower<MT>::value )
                   ?( ( IsUpper<MT>::value )?( sd.size() )
                                            :( ( IsStrictlyLower<MT>::value && sd.size() > 0UL )
                                               ?( sd.offset() + sd.size() - 1UL )
                                               :( sd.offset() + sd.size() ) ) )
                   :( ( IsUpper<MT>::value )?( left.columns() - column )
                                            :( left.columns() ) ) );

   return submatrix<AF>( left, sd.offset(), column, sd.size(), n ) * subvector<AF>( right, column, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/matrix multiplication.
// \ingroup submatrix
//
// \param vector The constant vector/matrix multiplication.
// \param args The runtime subvector arguments.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/matrix multiplication.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t... CSAs      // Compile time subvector arguments
        , typename VT         // Vector base type of the expression
        , typename... RSAs >  // Runtime subvector arguments
inline decltype(auto) subvector( const TVecMatMultExpr<VT>& vector, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using MT = RemoveReference_< RightOperand_< VectorType_<VT> > >;

   const SubvectorData<CSAs...> sd( args... );

   decltype(auto) left ( (~vector).leftOperand()  );
   decltype(auto) right( (~vector).rightOperand() );

   const size_t row( ( IsLower<MT>::value )
                     ?( ( !AF && IsStrictlyLower<MT>::value )?( sd.offset() + 1UL ):( sd.offset() ) )
                     :( 0UL ) );
   const size_t m( ( IsUpper<MT>::value )
                   ?( ( IsLower<MT>::value )?( sd.size() )
                                            :( ( IsStrictlyUpper<MT>::value && sd.size() > 0UL )
                                               ?( sd.offset() + sd.size() - 1UL )
                                               :( sd.offset() + sd.size() ) ) )
                   :( ( IsLower<MT>::value )?( right.rows() - row )
                                            :( right.rows() ) ) );

   return subvector<AF>( left, row, m ) * submatrix<AF>( right, row, sd.offset(), m, sd.size() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ROW)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given submatrix.
// \ingroup submatrix
//
// \param sm The submatrix containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the submatrix.
//
// This function returns an expression representing the specified row of the given submatrix.
*/
template< size_t I1           // Row index
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J            // Present index of the first column
        , size_t M            // Present number of rows
        , size_t N            // Present number of columns
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Submatrix<MT,AF,SO,DF,I2,J,M,N>& sm, RRAs... args )
{
   return subvector( row<I1+I2>( sm.operand(), args... ), J, N, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given constant submatrix.
// \ingroup submatrix
//
// \param sm The constant submatrix containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the submatrix.
//
// This function returns an expression representing the specified row of the given constant
// submatrix.
*/
template< size_t I1           // Row index
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J            // Present index of the first column
        , size_t M            // Present number of rows
        , size_t N            // Present number of columns
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const Submatrix<MT,AF,SO,DF,I2,J,M,N>& sm, RRAs... args )
{
   return subvector( row<I1+I2>( sm.operand(), args... ), J, N, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given temporary submatrix.
// \ingroup submatrix
//
// \param sm The temporary submatrix containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the submatrix.
//
// This function returns an expression representing the specified row of the given temporary
// submatrix.
*/
template< size_t I1           // Row index
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J            // Present index of the first column
        , size_t M            // Present number of rows
        , size_t N            // Present number of columns
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Submatrix<MT,AF,SO,DF,I2,J,M,N>&& sm, RRAs... args )
{
   return subvector( row<I1+I2>( sm.operand(), args... ), J, N, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given submatrix.
// \ingroup submatrix
//
// \param sm The submatrix containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the submatrix.
//
// This function returns an expression representing the specified row of the given submatrix.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Submatrix<MT,AF,SO,DF,CSAs...>& sm, RRAs... args )
{
   const RowData<CRAs...> rd( args... );
   const size_t index( rd.row() + sm.row() );

   return subvector( row( sm.operand(), index, args... ), sm.column(), sm.columns(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given constant submatrix.
// \ingroup submatrix
//
// \param sm The constant submatrix containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the submatrix.
//
// This function returns an expression representing the specified row of the given constant
// submatrix.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, RRAs... args )
{
   const RowData<CRAs...> rd( args... );
   const size_t index( rd.row() + sm.row() );

   return subvector( row( sm.operand(), index, args... ), sm.column(), sm.columns(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given temporary submatrix.
// \ingroup submatrix
//
// \param sm The temporary submatrix containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the submatrix.
//
// This function returns an expression representing the specified row of the given temporary
// submatrix.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Submatrix<MT,AF,SO,DF,CSAs...>&& sm, RRAs... args )
{
   const RowData<CRAs...> rd( args... );
   const size_t index( rd.row() + sm.row() );

   return subvector( row( sm.operand(), index, args... ), sm.column(), sm.columns(), unchecked );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (COLUMN)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given submatrix.
// \ingroup submatrix
//
// \param sm The submatrix containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the submatrix.
//
// This function returns an expression representing the specified column of the given submatrix.
*/
template< size_t I1           // Column index
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J            // Present index of the first column
        , size_t M            // Present number of rows
        , size_t N            // Present number of columns
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Submatrix<MT,AF,SO,DF,I2,J,M,N>& sm, RCAs... args )
{
   return subvector( column<I1+J>( sm.operand(), args... ), I2, M, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given constant submatrix.
// \ingroup submatrix
//
// \param sm The constant submatrix containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the submatrix.
//
// This function returns an expression representing the specified column of the given constant
// submatrix.
*/
template< size_t I1           // Column index
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J            // Present index of the first column
        , size_t M            // Present number of rows
        , size_t N            // Present number of columns
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const Submatrix<MT,AF,SO,DF,I2,J,M,N>& sm, RCAs... args )
{
   return subvector( column<I1+J>( sm.operand(), args... ), I2, M, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given temporary submatrix.
// \ingroup submatrix
//
// \param sm The temporary submatrix containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the submatrix.
//
// This function returns an expression representing the specified column of the given temporary
// submatrix.
*/
template< size_t I1           // Column index
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t I2           // Present index of the first row
        , size_t J            // Present index of the first column
        , size_t M            // Present number of rows
        , size_t N            // Present number of columns
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Submatrix<MT,AF,SO,DF,I2,J,M,N>&& sm, RCAs... args )
{
   return subvector( column<I1+J>( sm.operand(), args... ), I2, M, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given submatrix.
// \ingroup submatrix
//
// \param sm The submatrix containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the submatrix.
//
// This function returns an expression representing the specified column of the given submatrix.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Submatrix<MT,AF,SO,DF,CSAs...>& sm, RCAs... args )
{
   const ColumnData<CCAs...> cd( args... );
   const size_t index( cd.column() + sm.column() );

   return subvector( column( sm.operand(), index, args... ), sm.row(), sm.rows(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given constant submatrix.
// \ingroup submatrix
//
// \param sm The constant submatrix containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the submatrix.
//
// This function returns an expression representing the specified column of the given constant
// submatrix.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, RCAs... args )
{
   const ColumnData<CCAs...> cd( args... );
   const size_t index( cd.column() + sm.column() );

   return subvector( column( sm.operand(), index, args... ), sm.row(), sm.rows(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given temporary submatrix.
// \ingroup submatrix
//
// \param sm The temporary submatrix containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the submatrix.
//
// This function returns an expression representing the specified column of the given temporary
// submatrix.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the sparse submatrix
        , AlignmentFlag AF    // Present alignment flag
        , bool SO             // Storage order
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Submatrix<MT,AF,SO,DF,CSAs...>&& sm, RCAs... args )
{
   const ColumnData<CCAs...> cd( args... );
   const size_t index( cd.column() + sm.column() );

   return subvector( column( sm.operand(), index, args... ), sm.row(), sm.rows(), unchecked );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given submatrix.
// \ingroup submatrix
//
// \param sm The submatrix to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline void reset( Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   sm.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary submatrix.
// \ingroup submatrix
//
// \param sm The temporary submatrix to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline void reset( Submatrix<MT,AF,SO,DF,CSAs...>&& sm )
{
   sm.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row/column of the given submatrix.
// \ingroup submatrix
//
// \param sm The submatrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given submatrix to their
// default value. In case the given submatrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline void reset( Submatrix<MT,AF,SO,DF,CSAs...>& sm, size_t i )
{
   sm.reset( i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given matrix.
// \ingroup submatrix
//
// \param sm The matrix to be cleared.
// \return void
//
// Clearing a submatrix is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline void clear( Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   sm.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary matrix.
// \ingroup submatrix
//
// \param sm The temporary matrix to be cleared.
// \return void
//
// Clearing a submatrix is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline void clear( Submatrix<MT,AF,SO,DF,CSAs...>&& sm )
{
   sm.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense submatrix is in default state.
// \ingroup submatrix
//
// \param sm The dense submatrix to be tested for its default state.
// \return \a true in case the given dense submatrix is component-wise zero, \a false otherwise.
//
// This function checks whether the dense submatrix is in default state. For instance, in case
// the submatrix is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all submatrix elements are 0 and \a false in case any submatrix element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the dense matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isDefault( const Submatrix<MT,AF,SO,true,CSAs...>& sm )
{
   using blaze::isDefault;

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<(~sm).rows(); ++i )
         for( size_t j=0UL; j<(~sm).columns(); ++j )
            if( !isDefault<RF>( (~sm)(i,j) ) )
               return false;
   }
   else {
      for( size_t j=0UL; j<(~sm).columns(); ++j )
         for( size_t i=0UL; i<(~sm).rows(); ++i )
            if( !isDefault<RF>( (~sm)(i,j) ) )
               return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse submatrix is in default state.
// \ingroup submatrix
//
// \param sm The sparse submatrix to be tested for its default state.
// \return \a true in case the given sparse submatrix is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse submatrix is in default state. For instance, in case
// the submatrix is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all submatrix elements are 0 and \a false in case any submatrix element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the sparse matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isDefault( const Submatrix<MT,AF,SO,false,CSAs...>& sm )
{
   using blaze::isDefault;

   const size_t iend( ( SO == rowMajor)?( sm.rows() ):( sm.columns() ) );

   for( size_t i=0UL; i<iend; ++i ) {
      for( auto element=sm.cbegin(i); element!=sm.cend(i); ++element )
         if( !isDefault<RF>( element->value() ) ) return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given submatrix are intact.
// \ingroup submatrix
//
// \param sm The submatrix to be tested.
// \return \a true in case the given submatrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the submatrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( submatrix( A, 12UL, 13UL, 22UL, 33UL ) ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isIntact( const Submatrix<MT,AF,SO,DF,CSAs...>& sm ) noexcept
{
   return ( sm.row() + sm.rows() <= sm.operand().rows() &&
            sm.column() + sm.columns() <= sm.operand().columns() &&
            isIntact( sm.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is symmetric.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is symmetric, \a false if not.
//
// This function checks if the given submatrix is symmetric. The submatrix is considered to
// be symmetric if it is a square matrix whose transpose is equal to itself (\f$ A = A^T \f$). The
// following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isSymmetric( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isSymmetric( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsSymmetric<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isSymmetric( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is Hermitian.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is Hermitian, \a false if not.
//
// This function checks if the given submatrix is Hermitian. The submatrix is considered to
// be Hermitian if it is a square matrix whose transpose is equal to its conjugate transpose
// (\f$ A = \overline{A^T} \f$). The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isHermitian( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isHermitian( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsHermitian<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isHermitian( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is a lower triangular matrix.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is a lower triangular matrix, \a false if not.
//
// This function checks if the given submatrix is a lower triangular matrix. The matrix is
// considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        l_{1,0} & l_{1,1} & 0       & \cdots & 0       \\
                        l_{2,0} & l_{2,1} & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially lower triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isLower( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isLower( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsLower<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isLower( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is a lower unitriangular matrix.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is a lower unitriangular matrix, \a false if not.
//
// This function checks if the given submatrix is a lower unitriangular matrix. The matrix is
// considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 1       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 1       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 1      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isUniLower( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isUniLower( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsUniLower<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isUniLower( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is a strictly lower triangular matrix.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is a strictly lower triangular matrix, \a false if not.
//
// This function checks if the given submatrix is a strictly lower triangular matrix. The
// matrix is considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 0       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 0       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 0      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isStrictlyLower( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isStrictlyLower( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsStrictlyLower<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isStrictlyLower( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is an upper triangular matrix.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is an upper triangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is an upper triangular matrix. The matrix
// is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        u_{0,0} & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & u_{1,1} & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & u_{2,2} & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & u_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially upper triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isUpper( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsUpper<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isUpper( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is an upper unitriangular matrix.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is an upper unitriangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is an upper triangular matrix. The matrix
// is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 1       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 1       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 1       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isUniUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isUniUpper( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsUniUpper<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isUniUpper( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given submatrix is a strictly upper triangular matrix.
// \ingroup submatrix
//
// \param sm The submatrix to be checked.
// \return \a true if the submatrix is a strictly upper triangular matrix, \a false if not.
//
// This function checks if the given sparse submatrix is a strictly upper triangular matrix. The
// matrix is considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 0       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 0       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 0       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A( 32UL, 16UL );
   // ... Initialization

   auto sm = submatrix( A, 8UL, 8UL, 16UL, 16UL );

   if( isStrictlyUpper( sm ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isStrictlyUpper( const Submatrix<MT,AF,SO,DF,CSAs...>& sm )
{
   using BaseType = BaseType_< Submatrix<MT,AF,SO,DF,CSAs...> >;

   if( IsStrictlyUpper<MT>::value && sm.row() == sm.column() && sm.rows() == sm.columns() )
      return true;
   else return isStrictlyUpper( static_cast<const BaseType&>( sm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given matrix and submatrix represent the same observable state.
// \ingroup submatrix
//
// \param a The submatrix to be tested for its state.
// \param b The matrix to be tested for its state.
// \return \a true in case the submatrix and matrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given submatrix refers to the full given
// matrix and by that represents the same observable state. In this case, the function returns
// \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isSame( const Submatrix<MT,AF,SO,DF,CSAs...>& a, const Matrix<MT,SO>& b ) noexcept
{
   return ( isSame( a.operand(), ~b ) &&
            ( a.rows() == (~b).rows() ) &&
            ( a.columns() == (~b).columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given matrix and submatrix represent the same observable state.
// \ingroup submatrix
//
// \param a The matrix to be tested for its state.
// \param b The submatrix to be tested for its state.
// \return \a true in case the matrix and submatrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given submatrix refers to the full given
// matrix and by that represents the same observable state. In this case, the function returns
// \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , AlignmentFlag AF  // Alignment flag
        , bool DF           // Density flag
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isSame( const Matrix<MT,SO>& a, const Submatrix<MT,AF,SO,DF,CSAs...>& b ) noexcept
{
   return ( isSame( ~a, b.operand() ) &&
            ( (~a).rows() == b.rows() ) &&
            ( (~a).columns() == b.columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given submatrices represent the same observable state.
// \ingroup submatrix
//
// \param a The first submatrix to be tested for its state.
// \param b The second submatrix to be tested for its state.
// \return \a true in case the two submatrices share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given submatrices refer to exactly the
// same part of the same matrix. In case both submatrices represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the matrix of the left-hand side submatrix
        , AlignmentFlag AF1  // Alignment flag of the left-hand side submatrix
        , bool SO1           // Storage order of the left-hand side submatrix
        , bool DF1           // Density flag of the left-hand side submatrix
        , size_t... CSAs1    // Compile time submatrix arguments of the left-hand side submatrix
        , typename MT2       // Type of the matrix of the right-hand side submatrix
        , AlignmentFlag AF2  // Alignment flag of the right-hand side submatrix
        , bool SO2           // Storage order of the right-hand side submatrix
        , bool DF2           // Density flag of the right-hand side submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the right-hand side submatrix
inline bool isSame( const Submatrix<MT1,AF1,SO1,DF1,CSAs1...>& a,
                    const Submatrix<MT2,AF2,SO2,DF2,CSAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) &&
            ( a.row() == b.row() ) && ( a.column() == b.column() ) &&
            ( a.rows() == b.rows() ) && ( a.columns() == b.columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense submatrix.
// \ingroup submatrix
//
// \param sm The dense submatrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense submatrix by means of the specified matrix type or matrix
// inversion algorithm \c IF (see the InversionFlag documentation):

   \code
   invert<asLower>( A );     // Inversion of a lower triangular matrix
   invert<asUniUpper>( A );  // Inversion of an upper unitriangular matrix
   invert<byLU>( A );        // Inversion by means of an LU decomposition
   invert<byLLH>( A );       // Inversion by means of a Cholesky decomposition
   ...
   \endcode

// The matrix inversion fails if ...
//
//  - ... the given submatrix is not a square matrix;
//  - ... the given submatrix is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a sm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , size_t... CSAs >  // Compile time submatrix arguments
inline DisableIf_< HasMutableDataAccess<MT> > invert( Submatrix<MT,AF,SO,true,CSAs...>& sm )
{
   using RT = ResultType_< Submatrix<MT,AF,SO,true,CSAs...> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION  ( RT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( RT );

   RT tmp( sm );
   invert<IF>( tmp );
   sm = tmp;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a submatrix.
// \ingroup submatrix
//
// \param sm The target submatrix.
// \param i The row index of the element to be set.
// \param j The column index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename ET >     // Type of the element
inline bool trySet( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i <= sm.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j <= sm.columns(), "Invalid column access index" );

   return trySet( sm.operand(), sm.row()+i, sm.column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a submatrix.
// \ingroup submatrix
//
// \param sm The target submatrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename ET >     // Type of the element
inline bool tryAdd( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i <= sm.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j <= sm.columns(), "Invalid column access index" );

   return tryAdd( sm.operand(), sm.row()+i, sm.column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a submatrix.
// \ingroup submatrix
//
// \param sm The target submatrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename ET >     // Type of the element
inline bool trySub( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i <= sm.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j <= sm.columns(), "Invalid column access index" );

   return trySub( sm.operand(), sm.row()+i, sm.column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a submatrix.
// \ingroup submatrix
//
// \param sm The target submatrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename ET >     // Type of the element
inline bool tryMult( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i <= sm.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j <= sm.columns(), "Invalid column access index" );

   return tryMult( sm.operand(), sm.row()+i, sm.column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a submatrix.
// \ingroup submatrix
//
// \param sm The target submatrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename ET >     // Type of the element
inline bool tryDiv( const Submatrix<MT,AF,SO,DF,CSAs...>& sm, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i <= sm.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j <= sm.columns(), "Invalid column access index" );

   return tryDiv( sm.operand(), sm.row()+i, sm.column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                       const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to the band of a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be assigned.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                       const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryAssign( lhs.operand(), ~rhs, band + ptrdiff_t( lhs.column() - lhs.row() ),
                     lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a matrix to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryAssign( const Submatrix<MT1,AF,SO1,DF,CSAs...>& lhs,
                       const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return tryAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryAddAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to the band of
//        a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be added.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryAddAssign( lhs.operand(), ~rhs, band + ptrdiff_t( lhs.column() - lhs.row() ),
                        lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side matrix to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryAddAssign( const Submatrix<MT1,AF,SO1,DF,CSAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return tryAddAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool trySubAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return trySubAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to the band of
//        a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be subtracted.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool trySubAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return trySubAssign( lhs.operand(), ~rhs, band + ptrdiff_t( lhs.column() - lhs.row() ),
                        lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side matrix to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool trySubAssign( const Submatrix<MT1,AF,SO1,DF,CSAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return trySubAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be multiplied.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                           const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryMultAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to the band
//        of a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector to be multiplied.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                           const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryMultAssign( lhs.operand(), ~rhs, band + ptrdiff_t( lhs.column() - lhs.row() ),
                         lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the Schur product assignment of a matrix to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side matrix for the Schur product.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool trySchurAssign( const Submatrix<MT1,AF,SO1,DF,CSAs...>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   return trySchurAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector divisor.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryDivAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( TF || ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( !TF || ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryDivAssign( lhs.operand(), ~rhs, lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to the band of
//        a submatrix.
// \ingroup submatrix
//
// \param lhs The target left-hand side submatrix.
// \param rhs The right-hand side vector divisor.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryDivAssign( const Submatrix<MT,AF,SO,DF,CSAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.rows() - row ), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( ( (~rhs).size() <= lhs.columns() - column ), "Invalid number of columns" );

   return tryDivAssign( lhs.operand(), ~rhs, band + ptrdiff_t( lhs.column() - lhs.row() ),
                        lhs.row() + row, lhs.column() + column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given submatrix.
// \ingroup submatrix
//
// \param dm The submatrix to be derestricted.
// \return Submatrix without access restrictions.
//
// This function removes all restrictions on the data access to the given submatrix. It returns a
// submatrix that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t I          // Index of the first row
        , size_t J          // Index of the first column
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
inline decltype(auto) derestrict( Submatrix<MT,AF,SO,DF,I,J,M,N>& dm )
{
   return submatrix<AF,I,J,M,N>( derestrict( dm.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary submatrix.
// \ingroup submatrix
//
// \param dm The temporary submatrix to be derestricted.
// \return Submatrix without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary submatrix. It
// returns a submatrix that does provide the same interface but does not have any restrictions on
// the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF           // Density flag
        , size_t I          // Index of the first row
        , size_t J          // Index of the first column
        , size_t M          // Number of rows
        , size_t N >        // Number of columns
inline decltype(auto) derestrict( Submatrix<MT,AF,SO,DF,I,J,M,N>&& dm )
{
   return submatrix<AF,I,J,M,N>( derestrict( dm.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given submatrix.
// \ingroup submatrix
//
// \param dm The submatrix to be derestricted.
// \return Submatrix without access restrictions.
//
// This function removes all restrictions on the data access to the given submatrix. It returns a
// submatrix that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF >         // Density flag
inline decltype(auto) derestrict( Submatrix<MT,AF,SO,DF>& dm )
{
   return submatrix<AF>( derestrict( dm.operand() ), dm.row(), dm.column(), dm.rows(), dm.columns(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary submatrix.
// \ingroup submatrix
//
// \param dm The temporary submatrix to be derestricted.
// \return Submatrix without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary submatrix. It
// returns a submatrix that does provide the same interface but does not have any restrictions on
// the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO           // Storage order
        , bool DF >         // Density flag
inline decltype(auto) derestrict( Submatrix<MT,AF,SO,DF>&& dm )
{
   return submatrix<AF>( derestrict( dm.operand() ), dm.row(), dm.column(), dm.rows(), dm.columns(), unchecked );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct Size< Submatrix<MT,AF,SO,DF,I,J,M,N>, 0UL >
   : public PtrdiffT<M>
{};

template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct Size< Submatrix<MT,AF,SO,DF,I,J,M,N>, 1UL >
   : public PtrdiffT<N>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs >
struct IsRestricted< Submatrix<MT,AF,SO,DF,CSAs...> >
   : public IsRestricted<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, size_t... CSAs >
struct HasConstDataAccess< Submatrix<MT,AF,SO,true,CSAs...> >
   : public HasConstDataAccess<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, size_t... CSAs >
struct HasMutableDataAccess< Submatrix<MT,AF,SO,true,CSAs...> >
   : public HasMutableDataAccess<MT>
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
template< typename MT, bool SO, size_t... CSAs >
struct IsAligned< Submatrix<MT,aligned,SO,true,CSAs...> >
   : public TrueType
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsSymmetric< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsSymmetric<MT>::value && I == J && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsHermitian< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsHermitian<MT>::value && I == J && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsLower< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsLower<MT>::value && I == J && M == N ) ||
                          ( IsStrictlyLower<MT>::value && I == J+1UL && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsUniLower< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsUniLower<MT>::value && I == J && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsStrictlyLower< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsLower<MT>::value && I < J && M == N ) ||
                          ( IsStrictlyLower<MT>::value && I == J && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsUpper< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsUpper<MT>::value && I == J && M == N ) ||
                          ( IsStrictlyUpper<MT>::value && I+1UL == J && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsUniUpper< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsUniUpper<MT>::value && I == J && M == N ) >
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t I, size_t J, size_t M, size_t N >
struct IsStrictlyUpper< Submatrix<MT,AF,SO,DF,I,J,M,N> >
   : public BoolConstant< ( IsUpper<MT>::value && I > J && M == N ) ||
                          ( IsStrictlyUpper<MT>::value && I == J && M == N ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs1, size_t... CSAs2 >
struct SubmatrixTrait< Submatrix<MT,AF,SO,DF,CSAs1...>, CSAs2... >
{
   using Type = SubmatrixTrait_< ResultType_< Submatrix<MT,AF,SO,DF,CSAs1...> >, CSAs2... >;
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs, size_t... CRAs >
struct RowTrait< Submatrix<MT,AF,SO,DF,CSAs...>, CRAs... >
{
   using Type = RowTrait_< ResultType_< Submatrix<MT,AF,SO,DF,CSAs...> >, CRAs... >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs, size_t... CRAs >
struct RowsTrait< Submatrix<MT,AF,SO,DF,CSAs...>, CRAs... >
{
   using Type = RowsTrait_< ResultType_< Submatrix<MT,AF,SO,DF,CSAs...> >, CRAs... >;
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
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs, size_t... CCAs >
struct ColumnTrait< Submatrix<MT,AF,SO,DF,CSAs...>, CCAs... >
{
   using Type = ColumnTrait_< ResultType_< Submatrix<MT,AF,SO,DF,CSAs...> >, CCAs... >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs, size_t... CCAs >
struct ColumnsTrait< Submatrix<MT,AF,SO,DF,CSAs...>, CCAs... >
{
   using Type = ColumnsTrait_< ResultType_< Submatrix<MT,AF,SO,DF,CSAs...> >, CCAs... >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, AlignmentFlag AF, bool SO, bool DF, size_t... CSAs, ptrdiff_t... CBAs >
struct BandTrait< Submatrix<MT,AF,SO,DF,CSAs...>, CBAs... >
{
   using Type = BandTrait_< ResultType_< Submatrix<MT,AF,SO,DF,CSAs...> >, CBAs... >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
