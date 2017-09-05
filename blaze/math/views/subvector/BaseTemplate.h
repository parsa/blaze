//=================================================================================================
/*!
//  \file blaze/math/views/subvector/BaseTemplate.h
//  \brief Header file for the implementation of the Subvector base template
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_BASETEMPLATE_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_BASETEMPLATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsRowVector.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , bool TF          // Transpose flag
        , bool DF          // Density flag
        , size_t... SAs >  // Compile time subvector arguments
class SubvectorImpl
{};
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Reference to a specific subvector of a dense or sparse vector.
// \ingroup subvector
//
// The Subvector template represents a reference to a specific subvector of a dense or sparse
// vector primitive.
*/
template< typename VT          // Type of the vector
        , bool AF = unaligned  // Alignment flag
        , size_t... SAs >      // Compile time subvector arguments
using Subvector = SubvectorImpl< VT
                               , AF
                               , IsRowVector<VT>::value
                               , IsDenseVector<VT>::value
                               , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reference to a specific subvector of a dense vector.
// \ingroup subvector
//
// The DenseSubvector template represents a reference to a specific subvector of a dense vector
// primitive.
*/
template< typename VT          // Type of the vector
        , bool AF = unaligned  // Alignment flag
        , size_t... SAs >      // Compile time subvector arguments
using DenseSubvector = SubvectorImpl< VT
                                    , AF
                                    , IsRowVector<VT>::value
                                    , true
                                    , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reference to a specific subvector of a sparse vector.
// \ingroup subvector
//
// The SparseSubvector template represents a reference to a specific subvector of a sparse vector
// primitive.
*/
template< typename VT          // Type of the vector
        , bool AF = unaligned  // Alignment flag
        , size_t... SAs >      // Compile time subvector arguments
using SparseSubvector = SubvectorImpl< VT
                                     , AF
                                     , IsRowVector<VT>::value
                                     , false
                                     , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reference to a specific unaligned subvector of a vector.
// \ingroup subvector
//
// The UnalignedSubvector template represents a reference to a specific unaligned subvector of a
// vector primitive.
*/
template< typename VT      // Type of the vector
        , size_t... SAs >  // Compile time subvector arguments
using UnalignedSubvector = SubvectorImpl< VT
                                        , unaligned
                                        , IsRowVector<VT>::value
                                        , IsDenseVector<VT>::value
                                        , SAs... >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reference to a specific aligned subvector of a vector.
// \ingroup subvector
//
// The AlignedSubvector template represents a reference to a specific aligned subvector of a
// vector primitive.
*/
template< typename VT      // Type of the vector
        , size_t... SAs >  // Compile time subvector arguments
using AlignedSubvector = SubvectorImpl< VT
                                      , aligned
                                      , IsRowVector<VT>::value
                                      , IsDenseVector<VT>::value
                                      , SAs... >;
//*************************************************************************************************

} // namespace blaze

#endif
