//=================================================================================================
/*!
//  \file blaze/math/views/Band.h
//  \brief Header file for the implementation of the Band view
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

#ifndef _BLAZE_MATH_VIEWS_BAND_H_
#define _BLAZE_MATH_VIEWS_BAND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/DeclExpr.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSubmatrix.h>
#include <blaze/math/views/band/BaseTemplate.h>
#include <blaze/math/views/band/Dense.h>
#include <blaze/math/views/band/Sparse.h>
#include <blaze/math/views/Forward.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given matrix.
// \ingroup band
//
// \param matrix The matrix containing the band.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto ub1 = band<1L>( D );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto lb1 = band<-1L>( S );
   \endcode

// In case the band is not properly specified (i.e. if the specified index does not correspond
// to a valid band in the given matrix) a \a std::invalid_argument exception is thrown.
*/
template< ptrdiff_t I  // Band index
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Band<MT,I> band( Matrix<MT,SO>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return Band<MT,I>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given constant matrix.
// \ingroup band
//
// \param matrix The constant matrix containing the band.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto ub1 = band<1L>( D );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto lb1 = band<-1L>( S );
   \endcode

// In case the band is not properly specified (i.e. if the specified index does not correspond
// to a valid band in the given matrix) a \a std::invalid_argument exception is thrown.
*/
template< ptrdiff_t I  // Band index
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline const Band<const MT,I> band( const Matrix<MT,SO>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return Band<const MT,I>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given temporary matrix.
// \ingroup band
//
// \param matrix The temporary matrix containing the band.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given temporary
// matrix. In case the band is not properly specified (i.e. if the specified index does not
// correspond to a valid band in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< ptrdiff_t I  // Band index
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Band<MT,I> band( Matrix<MT,SO>&& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return Band<MT,I>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given matrix.
// \ingroup band
//
// \param matrix The matrix containing the band.
// \param index The band index.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto ub1 = band( D, 1L );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto lb1 = band( S, -1L );
   \endcode

// In case the band is not properly specified (i.e. if the specified index does not correspond
// to a valid band in the given matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Band<MT> band( Matrix<MT,SO>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return Band<MT>( ~matrix, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given constant matrix.
// \ingroup band
//
// \param matrix The constant matrix containing the band.
// \param index The band index.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto lb1 = band( D, 1L );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto ub1 = band( S, -1L );
   \endcode

// In case the band is not properly specified (i.e. if the specified index does not correspond
// to a valid band in the given matrix) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline const Band<const MT> band( const Matrix<MT,SO>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return Band<const MT>( ~matrix, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given temporary matrix.
// \ingroup band
//
// \param matrix The temporary matrix containing the band.
// \param index The band index.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given temporary
// matrix. In case the band is not properly specified (i.e. if the specified index does not
// correspond to a valid band in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Band<MT> band( Matrix<MT,SO>&& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return Band<MT>( ~matrix, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on the diagonal of the given matrix.
// \ingroup band
//
// \param matrix The matrix containing the diagonal.
// \return View on the diagonal of the matrix.
// \exception std::invalid_argument Invalid diagonal access index.
//
// This function returns an expression representing the diagonal of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the diagonal of the dense matrix D
   auto diag1 = diagonal( D );

   // Creating a view on the diagonal of the sparse matrix S
   auto diag2 = diagonal( S );
   \endcode

// In case the diagonal is not properly specified (i.e. in case the given matrix has zero rows
// or zero columns) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline decltype(auto) diagonal( Matrix<MT,SO>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<0L>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on the diagonal of the given constant matrix.
// \ingroup band
//
// \param matrix The constant matrix containing the diagonal.
// \return View on the diagonal of the matrix.
// \exception std::invalid_argument Invalid diagonal access index.
//
// This function returns an expression representing the diagonal of the given constant matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the diagonal of the dense matrix D
   auto diag1 = diagonal( D );

   // Creating a view on the diagonal of the sparse matrix S
   auto diag2 = diagonal( S );
   \endcode

// In case the diagonal is not properly specified (i.e. in case the given matrix has zero rows
// or zero columns) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline decltype(auto) diagonal( const Matrix<MT,SO>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<0L>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on the diagonal of the given temporary matrix.
// \ingroup band
//
// \param matrix The temporary matrix containing the diagonal.
// \return View on the diagonal of the matrix.
// \exception std::invalid_argument Invalid diagonal access index.
//
// This function returns an expression representing the diagonal of the given temporary matrix.
// In case the diagonal is not properly specified (i.e. in case the given matrix has zero rows
// or zero columns) a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline decltype(auto) diagonal( Matrix<MT,SO>&& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<0L>( ~matrix );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given band.
// \ingroup band
//
// \param b The band containing the subvector.
// \return View on the specified subvector of the band.
//
// This function returns an expression representing the specified subvector of the given band.
*/
template< bool AF         // Alignment flag
        , size_t I1       // Index of the first subvector element
        , size_t N        // Size of the subvector
        , typename MT     // Type of the matrix
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , bool MF         // Multiplication flag
        , ptrdiff_t I2 >  // Band index
inline decltype(auto) subvector( BandImpl<MT,TF,DF,MF,I2>& b )
{
   BLAZE_FUNCTION_TRACE;

   constexpr size_t row   ( I2 >= 0L ? 0UL : -I2 );
   constexpr size_t column( I2 >= 0L ?  I2 : 0UL );

   return diagonal( submatrix<AF,row+I1,column+I1,N,N>( b.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given band.
// \ingroup band
//
// \param b The band containing the subvector.
// \return View on the specified subvector of the band.
//
// This function returns an expression representing the specified subvector of the given band.
*/
template< bool AF      // Alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) subvector( BandImpl<MT,TF,DF,MF>& b )
{
   BLAZE_FUNCTION_TRACE;

   return diagonal( submatrix<AF>( b.operand(), b.row()+I, b.column()+I, N, N ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given band.
// \ingroup band
//
// \param b The band containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the band.
//
// This function returns an expression representing the specified subvector of the given band.
*/
template< bool AF             // Alignment flag
        , typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline decltype(auto) subvector( BandImpl<MT,TF,DF,MF,BAs...>& b, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return diagonal( submatrix<AF>( b.operand(), b.row()+index, b.column()+index, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given constant band.
// \ingroup band
//
// \param b The constant band containing the subvector.
// \return View on the specified subvector of the constant band.
//
// This function returns an expression representing the specified subvector of the given constant
// band.
*/
template< bool AF         // Alignment flag
        , size_t I1       // Index of the first subvector element
        , size_t N        // Size of the subvector
        , typename MT     // Type of the matrix
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , bool MF         // Multiplication flag
        , ptrdiff_t I2 >  // Band index
inline decltype(auto) subvector( const BandImpl<MT,TF,DF,MF,I2>& b )
{
   BLAZE_FUNCTION_TRACE;

   constexpr size_t row   ( I2 >= 0L ? 0UL : -I2 );
   constexpr size_t column( I2 >= 0L ?  I2 : 0UL );

   return diagonal( submatrix<AF,row+I1,column+I1,N,N>( b.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given constant band.
// \ingroup band
//
// \param b The constant band containing the subvector.
// \return View on the specified subvector of the constant band.
//
// This function returns an expression representing the specified subvector of the given constant
// band.
*/
template< bool AF      // Alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) subvector( const BandImpl<MT,TF,DF,MF>& b )
{
   BLAZE_FUNCTION_TRACE;

   return diagonal( submatrix<AF>( b.operand(), b.row()+I, b.column()+I, N, N ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given constant band.
// \ingroup band
//
// \param b The constant band containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the constant band.
//
// This function returns an expression representing the specified subvector of the given constant
// band.
*/
template< bool AF             // Alignment flag
        , typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline decltype(auto) subvector( const BandImpl<MT,TF,DF,MF,BAs...>& b, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return diagonal( submatrix<AF>( b.operand(), b.row()+index, b.column()+index, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given temporary band.
// \ingroup band
//
// \param b The temporary band containing the subvector.
// \return View on the specified subvector of the temporary band.
//
// This function returns an expression representing the specified subvector of the given temporary
// band.
*/
template< bool AF         // Alignment flag
        , size_t I1       // Index of the first subvector element
        , size_t N        // Size of the subvector
        , typename MT     // Type of the matrix
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , bool MF         // Multiplication flag
        , ptrdiff_t I2 >  // Band index
inline decltype(auto) subvector( BandImpl<MT,TF,DF,MF,I2>&& b )
{
   BLAZE_FUNCTION_TRACE;

   constexpr size_t row   ( I2 >= 0L ? 0UL : -I2 );
   constexpr size_t column( I2 >= 0L ?  I2 : 0UL );

   return diagonal( submatrix<AF,row+I1,column+I1,N,N>( b.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given temporary band.
// \ingroup band
//
// \param b The temporary band containing the subvector.
// \return View on the specified subvector of the temporary band.
//
// This function returns an expression representing the specified subvector of the given temporary
// band.
*/
template< bool AF      // Alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) subvector( BandImpl<MT,TF,DF,MF>&& b )
{
   BLAZE_FUNCTION_TRACE;

   return diagonal( submatrix<AF>( b.operand(), b.row()+I, b.column()+I, N, N ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given temporary band.
// \ingroup band
//
// \param b The temporary band containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the temporary band.
//
// This function returns an expression representing the specified subvector of the given temporary
// band.
*/
template< bool AF             // Alignment flag
        , typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline decltype(auto) subvector( BandImpl<MT,TF,DF,MF,BAs...>&& b, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return diagonal( submatrix<AF>( b.operand(), b.row()+index, b.column()+index, size, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/matrix addition.
// \ingroup band
//
// \param matrix The constant matrix/matrix addition.
// \return View on the specified band of the addition.
//
// This function returns an expression representing the specified band of the given matrix/matrix
// addition.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMatAddExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<I>( (~matrix).leftOperand() ) + band<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/matrix addition.
// \ingroup band
//
// \param matrix The constant matrix/matrix addition.
// \param index The band index.
// \return View on the specified band of the addition.
//
// This function returns an expression representing the specified band of the given matrix/matrix
// addition.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMatAddExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).leftOperand(), index ) + band( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/matrix subtraction.
// \ingroup band
//
// \param matrix The constant matrix/matrix subtraction.
// \return View on the specified band of the subtraction.
//
// This function returns an expression representing the specified band of the given matrix/matrix
// subtraction.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMatSubExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<I>( (~matrix).leftOperand() ) - band<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/matrix subtraction.
// \ingroup band
//
// \param matrix The constant matrix/matrix subtraction.
// \param index The band index.
// \return View on the specified band of the subtraction.
//
// This function returns an expression representing the specified band of the given matrix/matrix
// subtraction.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMatSubExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).leftOperand(), index ) - band( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given Schur product.
// \ingroup band
//
// \param matrix The constant Schur product.
// \return View on the specified band of the Schur product.
//
// This function returns an expression representing the specified band of the given Schur product.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const SchurExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<I>( (~matrix).leftOperand() ) * band<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given Schur product.
// \ingroup band
//
// \param matrix The constant Schur product.
// \param index The band index.
// \return View on the specified band of the Schur product.
//
// This function returns an expression representing the specified band of the given Schur product.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const SchurExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).leftOperand(), index ) * band( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given outer product.
// \ingroup band
//
// \param matrix The constant outer product.
// \return View on the specified band of the outer product.
//
// This function returns an expression representing the specified band of the given outer product.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const VecTVecMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   decltype(auto) leftOperand ( (~matrix).leftOperand()  );
   decltype(auto) rightOperand( (~matrix).rightOperand() );

   const size_t row   ( I <  0L ? -I : 0UL );
   const size_t column( I >= 0L ?  I : 0UL );
   const size_t size  ( min( leftOperand.size() - row, rightOperand.size() - column ) );

   return transTo<defaultTransposeFlag>( subvector( leftOperand , row   , size ) ) *
          transTo<defaultTransposeFlag>( subvector( rightOperand, column, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given outer product.
// \ingroup band
//
// \param matrix The constant outer product.
// \param index The band index.
// \return View on the specified band of the outer product.
//
// This function returns an expression representing the specified band of the given outer product.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const VecTVecMultExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   decltype(auto) leftOperand ( (~matrix).leftOperand()  );
   decltype(auto) rightOperand( (~matrix).rightOperand() );

   const size_t row   ( index <  0L ? -index : 0UL );
   const size_t column( index >= 0L ?  index : 0UL );
   const size_t size  ( min( leftOperand.size() - row, rightOperand.size() - column ) );

   return transTo<defaultTransposeFlag>( subvector( leftOperand , row   , size ) ) *
          transTo<defaultTransposeFlag>( subvector( rightOperand, column, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/scalar multiplication.
// \ingroup band
//
// \param matrix The constant matrix/scalar multiplication.
// \return View on the specified band of the multiplication.
//
// This function returns an expression representing the specified band of the given matrix/scalar
// multiplication.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<I>( (~matrix).leftOperand() ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/scalar multiplication.
// \ingroup band
//
// \param matrix The constant matrix/scalar multiplication.
// \param index The band index.
// \return View on the specified band of the multiplication.
//
// This function returns an expression representing the specified band of the given matrix/scalar
// multiplication.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatScalarMultExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).leftOperand(), index ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/scalar division.
// \ingroup band
//
// \param matrix The constant matrix/scalar division.
// \return View on the specified band of the division.
//
// This function returns an expression representing the specified band of the given matrix/scalar
// division.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatScalarDivExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<I>( (~matrix).leftOperand() ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/scalar division.
// \ingroup band
//
// \param matrix The constant matrix/scalar division.
// \param index The band index.
// \return View on the specified band of the division.
//
// This function returns an expression representing the specified band of the given matrix/scalar
// division.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatScalarDivExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).leftOperand(), index ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given unary matrix map operation.
// \ingroup band
//
// \param matrix The constant unary matrix map operation.
// \return View on the specified band of the unary map operation.
//
// This function returns an expression representing the specified band of the given unary matrix
// map operation.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( band<I>( (~matrix).operand() ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given unary matrix map operation.
// \ingroup band
//
// \param matrix The constant unary matrix map operation.
// \param index The band index.
// \return View on the specified band of the unary map operation.
//
// This function returns an expression representing the specified band of the given unary matrix
// map operation.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMapExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return map( band( (~matrix).operand(), index ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given binary matrix map operation.
// \ingroup band
//
// \param matrix The constant binary matrix map operation.
// \return View on the specified band of the binary map operation.
//
// This function returns an expression representing the specified band of the given binary matrix
// map operation.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMatMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( band<I>( (~matrix).leftOperand() ),
               band<I>( (~matrix).rightOperand() ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given binary matrix map operation.
// \ingroup band
//
// \param matrix The constant binary matrix map operation.
// \param index The band index.
// \return View on the specified band of the binary map operation.
//
// This function returns an expression representing the specified band of the given binary matrix
// map operation.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatMatMapExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return map( band( (~matrix).leftOperand(), index ),
               band( (~matrix).rightOperand(), index ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix evaluation operation.
// \ingroup band
//
// \param matrix The constant matrix evaluation operation.
// \return View on the specified band of the evaluation operation.
//
// This function returns an expression representing the specified band of the given matrix
// evaluation operation.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatEvalExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return eval( band<I>( (~matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix evaluation operation.
// \ingroup band
//
// \param matrix The constant matrix evaluation operation.
// \param index The band index.
// \return View on the specified band of the evaluation operation.
//
// This function returns an expression representing the specified band of the given matrix
// evaluation operation.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatEvalExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return eval( band( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix serialization operation.
// \ingroup band
//
// \param matrix The constant matrix serialization operation.
// \return View on the specified band of the serialization operation.
//
// This function returns an expression representing the specified band of the given matrix
// serialization operation.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatSerialExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return serial( band<I>( (~matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix serialization operation.
// \ingroup band
//
// \param matrix The constant matrix serialization operation.
// \param index The band index.
// \return View on the specified band of the serialization operation.
//
// This function returns an expression representing the specified band of the given matrix
// serialization operation.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatSerialExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return serial( band( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix declaration operation.
// \ingroup band
//
// \param matrix The constant matrix declaration operation.
// \return View on the specified band of the declaration operation.
//
// This function returns an expression representing the specified band of the given matrix
// declaration operation.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const DeclExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<I>( (~matrix).operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix declaration operation.
// \ingroup band
//
// \param matrix The constant matrix declaration operation.
// \param index The band index.
// \return View on the specified band of the declaration operation.
//
// This function returns an expression representing the specified band of the given matrix
// declaration operation.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const DeclExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).operand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix transpose operation.
// \ingroup band
//
// \param matrix The constant matrix transpose operation.
// \return View on the specified band of the transpose operation.
//
// This function returns an expression representing the specified band of the given matrix
// transpose operation.
*/
template< ptrdiff_t I    // Band index
        , typename MT >  // Type of the matrix
inline decltype(auto) band( const MatTransExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return band<-I>( (~matrix).operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix transpose operation.
// \ingroup band
//
// \param matrix The constant matrix declaration operation.
// \param index The band index.
// \return View on the specified band of the declaration operation.
//
// This function returns an expression representing the specified band of the given matrix
// declaration operation.
*/
template< typename MT >  // Type of the matrix
inline decltype(auto) band( const MatTransExpr<MT>& matrix, ptrdiff_t index )
{
   BLAZE_FUNCTION_TRACE;

   return band( (~matrix).operand(), -index );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BAND OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given band.
// \ingroup band
//
// \param band The band to be resetted.
// \return void
*/
template< typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline void reset( BandImpl<MT,TF,DF,MF,BAs...>& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary band.
// \ingroup band
//
// \param band The temporary band to be resetted.
// \return void
*/
template< typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline void reset( BandImpl<MT,TF,DF,MF,BAs...>&& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given band.
// \ingroup band
//
// \param band The band to be cleared.
// \return void
//
// Clearing a band is equivalent to resetting it via the reset() function.
*/
template< typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline void clear( BandImpl<MT,TF,DF,MF,BAs...>& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary band.
// \ingroup band
//
// \param band The temporary band to be cleared.
// \return void
//
// Clearing a band is equivalent to resetting it via the reset() function.
*/
template< typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline void clear( BandImpl<MT,TF,DF,MF,BAs...>&& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense band is in default state.
// \ingroup band
//
// \param band The dense band to be tested for its default state.
// \return \a true in case the given dense band is component-wise zero, \a false otherwise.
//
// This function checks whether the dense band is in default state. For instance, in case the
// band is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all band elements are 0 and \a false in case any band element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( band( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( band( A, 0UL ) ) ) { ... }
   \endcode
*/
template< bool RF             // Relaxation flag
        , typename MT         // Type of the dense matrix
        , bool TF             // Transpose flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline bool isDefault( const BandImpl<MT,TF,true,MF,BAs...>& band )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<band.size(); ++i )
      if( !isDefault<RF>( band[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse band is in default state.
// \ingroup band
//
// \param band The sparse band to be tested for its default state.
// \return \a true in case the given band is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse band is in default state. For instance, in case the
// band is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all band elements are 0 and \a false in case any band element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( band( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( band( A, 0UL ) ) ) { ... }
   \endcode
*/
template< bool RF             // Relaxation flag
        , typename MT         // Type of the sparse matrix
        , bool TF             // Transpose flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline bool isDefault( const BandImpl<MT,TF,false,MF,BAs...>& band )
{
   using blaze::isDefault;

   for( const auto& element : band )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given band are intact.
// \ingroup band
//
// \param band The band to be tested.
// \return \a true in case the given band's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the band are intact, i.e. if its state is valid.
// In case the invariants are intact, the function returns \a true, else it will return \a false.
// The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( band( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , bool MF             // Multiplication flag
        , ptrdiff_t... BAs >  // Compile time band arguments
inline bool isIntact( const BandImpl<MT,TF,DF,MF,BAs...>& band ) noexcept
{
   const ptrdiff_t index( band.band() );

   return ( ( index >= 0L || size_t( -index ) < band.operand().rows()    ) &&
            ( index <  0L || size_t(  index ) < band.operand().columns() ) &&
            isIntact( band.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two regular bands.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of two regular
// bands. In case both bands represent the same observable state, the function returns \a true,
// otherwise it returns \a false.
*/
template< typename MT1         // Type of the matrix of the left-hand side band
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... BAs1    // Compile time band arguments of the left-hand side band
        , typename MT2         // Type of the matrix of the right-hand side band
        , ptrdiff_t... BAs2 >  // Compile time band arguments of the right-hand side band
inline DisableIf_< Or< IsSubmatrix<MT1>, IsSubmatrix<MT2> >, bool >
   isSame_backend( const BandImpl<MT1,TF,DF,MF,BAs1...>& a,
                   const BandImpl<MT2,TF,DF,MF,BAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) && ( a.band() == b.band() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the left band being a band on a submatrix.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the left
// band being a band on a submatrix. In case both bands represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1         // Type of the submatrix of the left-hand side band
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... BAs1    // Compile time band arguments of the left-hand side band
        , typename MT2         // Type of the matrix of the right-hand side band
        , ptrdiff_t... BAs2 >  // Compile time band arguments of the right-hand side band
inline EnableIf_< And< IsSubmatrix<MT1>, Not< IsSubmatrix<MT2> > >, bool >
   isSame_backend( const BandImpl<MT1,TF,DF,MF,BAs1...>& a,
                   const BandImpl<MT2,TF,DF,MF,BAs2...>& b ) noexcept
{
   return ( isSame( a.operand().operand(), b.operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() + a.operand().row() == b.row() ) &&
            ( a.column() + a.operand().column() == b.column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the right band being a band on a submatrix.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the right
// band being a band on a submatrix. In case both bands represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1         // Type of the matrix of the left-hand side band
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... BAs1    // Compile time band arguments of the left-hand side band
        , typename MT2         // Type of the submatrix of the right-hand side band
        , ptrdiff_t... BAs2 >  // Compile time band arguments of the right-hand side band
inline EnableIf_< And< Not< IsSubmatrix<MT1> >, IsSubmatrix<MT2> >, bool >
   isSame_backend( const BandImpl<MT1,TF,DF,MF,BAs1...>& a,
                   const BandImpl<MT2,TF,DF,MF,BAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() == b.row() + b.operand().row() ) &&
            ( a.column() == b.column() + b.operand().column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two bands on submatrices.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of both bands
// being bands on submatrices. In case both bands represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1         // Type of the submatrix of the left-hand side band
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... BAs1    // Compile time band arguments of the left-hand side band
        , typename MT2         // Type of the submatrix of the right-hand side band
        , ptrdiff_t... BAs2 >  // Compile time band arguments of the right-hand side band
inline EnableIf_< And< IsSubmatrix<MT1>, IsSubmatrix<MT2> >, bool >
   isSame_backend( const BandImpl<MT1,TF,DF,MF,BAs1...>& a,
                   const BandImpl<MT2,TF,DF,MF,BAs2...>& b ) noexcept
{
   return ( isSame( a.operand().operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() + a.operand().row() == b.row() + b.operand().row() ) &&
            ( a.column() + a.operand().column() == b.column() + b.operand().column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given bands represent the same observable state.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given bands refer to exactly the same
// range of the same matrix. In case both bands represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1         // Type of the matrix of the left-hand side band
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... BAs1    // Compile time band arguments of the left-hand side band
        , typename MT2         // Type of the matrix of the right-hand side band
        , ptrdiff_t... BAs2 >  // Compile time band arguments of the right-hand side band
inline bool isSame( const BandImpl<MT1,TF,DF,MF,BAs1...>& a,
                    const BandImpl<MT2,TF,DF,MF,BAs2...>& b ) noexcept
{
   return isSame_backend( a, b );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , bool MF           // Multiplication flag
        , ptrdiff_t... BAs  // Compile time band arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryAssign( const BandImpl<MT,TF,DF,MF,BAs...>& lhs,
                       const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.operand(), ~rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , bool MF           // Multiplication flag
        , ptrdiff_t... BAs  // Compile time band arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryAddAssign( const BandImpl<MT,TF,DF,MF,BAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.operand(), ~rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , bool MF           // Multiplication flag
        , ptrdiff_t... BAs  // Compile time band arguments
        , typename VT >     // Type of the right-hand side vector
inline bool trySubAssign( const BandImpl<MT,TF,DF,MF,BAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.operand(), ~rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , bool MF           // Multiplication flag
        , ptrdiff_t... BAs  // Compile time band arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryMultAssign( const BandImpl<MT,TF,DF,MF,BAs...>& lhs,
                           const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.operand(), ~rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , bool MF           // Multiplication flag
        , ptrdiff_t... BAs  // Compile time band arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryDivAssign( const BandImpl<MT,TF,DF,MF,BAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryDivAssign( lhs.operand(), ~rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given band.
// \ingroup band
//
// \param b The band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given band. It returns a
// band object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT    // Type of the matrix
        , bool TF        // Transpose flag
        , bool DF        // Density flag
        , bool MF        // Multiplication flag
        , ptrdiff_t I >  // Band index
inline decltype(auto) derestrict( BandImpl<MT,TF,DF,MF,I>& b )
{
   return band<I>( derestrict( b.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary band.
// \ingroup band
//
// \param b The temporary band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary band. It
// returns a band object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT    // Type of the matrix
        , bool TF        // Transpose flag
        , bool DF        // Density flag
        , bool MF        // Multiplication flag
        , ptrdiff_t I >  // Band index
inline decltype(auto) derestrict( BandImpl<MT,TF,DF,MF,I>&& b )
{
   return band<I>( derestrict( b.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given band.
// \ingroup band
//
// \param b The band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given band. It returns a
// band object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) derestrict( BandImpl<MT,TF,DF,MF>& b )
{
   return band( derestrict( b.operand() ), b.band() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary band.
// \ingroup band
//
// \param b The temporary band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary band. It
// returns a band object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) derestrict( BandImpl<MT,TF,DF,MF>&& b )
{
   return band( derestrict( b.operand() ), b.band() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t... BAs >
struct IsRestricted< BandImpl<MT,TF,DF,MF,BAs...> >
   : public BoolConstant< IsRestricted<MT>::value >
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
template< typename MT, bool TF, bool MF, ptrdiff_t... BAs >
struct HasConstDataAccess< BandImpl<MT,TF,true,MF,BAs...> >
   : public BoolConstant< HasConstDataAccess<MT>::value >
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
template< typename MT, bool TF, bool MF, ptrdiff_t... BAs >
struct HasMutableDataAccess< BandImpl<MT,TF,true,MF,BAs...> >
   : public BoolConstant< HasMutableDataAccess<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISOPPOSEDVIEW SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t... BAs >
struct IsOpposedView< BandImpl<MT,TF,DF,MF,BAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t... BAs, size_t... CSAs >
struct SubvectorTrait< BandImpl<MT,TF,DF,MF,BAs...>, CSAs... >
{
   using Type = SubvectorTrait_< ResultType_< BandImpl<MT,TF,DF,MF,BAs...> >, CSAs... >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
