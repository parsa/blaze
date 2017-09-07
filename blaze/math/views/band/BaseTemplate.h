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
        , ptrdiff_t... BAs >  // Compile time band arguments
class BandImpl
{};
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief View on a specific band of a dense or sparse matrix.
// \ingroup band
//
// The Band template represents a reference to a specific band of a dense or sparse matrix primitive.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BAs >  // Compile time band arguments
using Band = BandImpl< MT
                     , defaultTransposeFlag
                     , IsDenseMatrix<MT>::value
                     , IsMatMatMultExpr<MT>::value
                     , BAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific band of a dense matrix.
// \ingroup band
//
// The DenseBand template represents a reference to a specific band of a dense matrix primitive.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BAs >  // Compile time band arguments
using DenseBand = BandImpl< MT
                          , defaultTransposeFlag
                          , true
                          , IsMatMatMultExpr<MT>::value
                          , BAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific band of a sparse matrix.
// \ingroup band
//
// The SparseBand template represents a reference to a specific band of a sparse matrix primitive.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BAs >  // Compile time band arguments
using SparseBand = BandImpl< MT
                           , defaultTransposeFlag
                           , false
                           , IsMatMatMultExpr<MT>::value
                           , BAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific diagonal of a dense or sparse matrix.
// \ingroup diagonal
//
// The DenseDiagonal template represents a reference to a specific diagonal of a dense or sparse
// matrix primitive.
*/
template< typename MT >
using Diagonal = Band<MT,0L>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific diagonal of a dense matrix.
// \ingroup diagonal
//
// The DenseDiagonal template represents a reference to a specific diagonal of a dense matrix
// primitive.
*/
template< typename MT >   // Type of the matrix
using DenseDiagonal = DenseBand<MT,0L>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific diagonal of a sparse matrix.
// \ingroup diagonal
//
// The SparseDiagonal template represents a reference to a specific diagonal of a sparse matrix
// primitive.
*/
template< typename MT >  // Type of the matrix
using SparseDiagonal = SparseBand<MT,0L>;
//*************************************************************************************************

} // namespace blaze

#endif
