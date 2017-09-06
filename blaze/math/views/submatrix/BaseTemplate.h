//=================================================================================================
/*!
//  \file blaze/math/views/submatrix/BaseTemplate.h
//  \brief Header file for the implementation of the Submatrix base template
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

#ifndef _BLAZE_MATH_VIEWS_SUBMATRIX_BASETEMPLATE_H_
#define _BLAZE_MATH_VIEWS_SUBMATRIX_BASETEMPLATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
template< typename MT      // Type of the matrix
        , bool AF          // Alignment flag
        , bool SO          // Storage order
        , bool DF          // Density flag
        , size_t... SAs >  // Compile time submatrix arguments
class SubmatrixImpl
{};
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief View on a specific submatrix of a dense or sparse matrix.
// \ingroup submatrix
//
// The Submatrix template represents a reference to a specific submatrix of a dense or sparse
// matrix primitive.
*/
template< typename MT          // Type of the matrix
        , bool AF = unaligned  // Alignment flag
        , size_t... SAs >      // Compile time submatrix arguments
using Submatrix = SubmatrixImpl< MT
                               , AF
                               , IsColumnMajorMatrix<MT>::value
                               , IsDenseMatrix<MT>::value
                               , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific submatrix of a dense matrix.
// \ingroup submatrix
//
// The DenseSubmatrix template represents a reference to a specific submatrix of a dense matrix
// primitive.
*/
template< typename MT          // Type of the matrix
        , bool AF = unaligned  // Alignment flag
        , size_t... SAs >      // Compile time submatrix arguments
using DenseSubmatrix = SubmatrixImpl< MT
                                    , AF
                                    , IsColumnMajorMatrix<MT>::value
                                    , true
                                    , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific submatrix of a sparse matrix.
// \ingroup submatrix
//
// The SparseSubmatrix template represents a reference to a specific submatrix of a sparse matrix
// primitive.
*/
template< typename MT          // Type of the matrix
        , bool AF = unaligned  // Alignment flag
        , size_t... SAs >      // Compile time submatrix arguments
using SparseSubmatrix = SubmatrixImpl< MT
                                     , AF
                                     , IsColumnMajorMatrix<MT>::value
                                     , true
                                     , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific unaligned submatrix of a matrix.
// \ingroup submatrix
//
// The UnalignedSubmatrix template represents a reference to a specific unaligned submatrix of a
// matrix primitive.
*/
template< typename MT      // Type of the matrix
        , size_t... SAs >  // Compile time submatrix arguments
using UnalignedSubmatrix = SubmatrixImpl< MT
                                        , unaligned
                                        , IsColumnMajorMatrix<MT>::value
                                        , IsDenseMatrix<MT>::value
                                        , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief View on a specific aligned submatrix of a matrix.
// \ingroup submatrix
//
// The AlignedSubmatrix template represents a reference to a specific aligned submatrix of a
// matrix primitive.
*/
template< typename MT      // Type of the matrix
        , size_t... SAs >  // Compile time submatrix arguments
using AlignedSubmatrix = SubmatrixImpl< MT
                                      , aligned
                                      , IsColumnMajorMatrix<MT>::value
                                      , IsDenseMatrix<MT>::value
                                      , SAs... >;
//*************************************************************************************************

} // namespace blaze

#endif
