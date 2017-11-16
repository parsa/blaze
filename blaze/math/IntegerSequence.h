//=================================================================================================
/*!
//  \file blaze/math/IndexSequence.h
//  \brief Header file for the std::index_sequence aliases
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

#ifndef _BLAZE_MATH_INDEXSEQUENCE_H_
#define _BLAZE_MATH_INDEXSEQUENCE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::integer_sequence
// \brief Integer sequence type of the Blaze library.
// \ingroup math
*/
using std::integer_sequence;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::index_sequence
// \brief Index sequence type of the Blaze library.
// \ingroup math
*/
using std::index_sequence;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::make_integer_sequence
// \brief Import of the std::make_integer_sequence alias template into the Blaze namespace.
// \ingroup math
*/
using std::make_integer_sequence;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::make_index_sequence
// \brief Import of the std::make_index_sequence alias template into the Blaze namespace.
// \ingroup math
*/
using std::make_index_sequence;
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Shifts the given index sequence by a given offset.
// \ingroup math
//
// \param sequence The given index sequence
// \return The shifted index sequence.
*/
template< size_t Offset   // The offset for the shift operation
        , size_t... Is >  // The sequence of indices
constexpr decltype(auto) shift( std::index_sequence<Is...> sequence )
{
   UNUSED_PARAMETER( sequence );

   return std::index_sequence< ( Is + Offset )... >();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the setup of shifted index sequences.
// \ingroup math
//
// The make_shifted_index_sequence alias template provides a convenient way to create index
// sequences with specific initial index and a specific number of indices. The following code
// example demonstrates the use of make_shifted_index_sequence:

   \code
   // Creating the index sequence <2,3,4,5,6>
   using Type = make_shifted_index_sequence<2UL,5UL>;
   \endcode
*/
template< size_t Offset  // The offset of the index sequence
        , size_t N >     // The total number of indices in the index sequence
using make_shifted_index_sequence = decltype( shift<Offset>( make_index_sequence<N>() ) );
//*************************************************************************************************

} // namespace blaze

#endif
