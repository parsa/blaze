//=================================================================================================
/*!
//  \file blaze/math/views/Check.h
//  \brief Header file for the blaze::checked and blaze::unchecked instances
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_CHECK_H_
#define _BLAZE_MATH_VIEWS_CHECK_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Template for the blaze::checked and blaze::unchecked instances.
// \ingroup math
//
// blaze::Check is the template for the blaze::checked and blaze::unchecked instance, which is
// an optional token for the creation of views. It can be used to enforce or skip all runtime
// checks during the creation of a view (subvectors, submatrices, rows, columns, bands, ...).
*/
template< bool C >
struct Check
{
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   constexpr Check() = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TYPE ALIASES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type of the blaze::checked instance.
// \ingroup math
//
// blaze::Checked is the type of the blaze::checked instance, which is an optional token for the
// creation of views. It can be used to enforce runtime checks during the creation of a view
// (subvectors, submatrices, rows, columns, bands, ...).
*/
using Checked = Check<true>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Type of the blaze::unchecked instance.
// \ingroup math
//
// blaze::Unchecked is the type of the blaze::unchecked instance, which is an optional token for
// the creation of views. It can be used to skip all runtime checks during the creation of a view
// (subvectors, submatrices, rows, columns, bands, ...).
*/
using Unchecked = Check<false>;
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL CHECK INSTANCES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global Checked instance.
// \ingroup math
//
// The blaze::checked instance is an optional token for the creation of views. It can be used
// used to enforce runtime checks during the creation of a view (subvectors, submatrices, rows,
// columns, bands, ...). The following example demonstrates the setup of a checked subvector:

   \code
   blaze::DynamicVector<int> v( 100UL );
   auto sv = subvector( v, 10UL, 20UL, checked );  // Creating an checked subvector
   \endcode
*/
constexpr Checked checked;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global Unchecked instance.
// \ingroup math
//
// The blaze::unchecked instance is an optional token for the creation of views. It can be
// used to skip all runtime checks during the creation of a view (subvectors, submatrices, rows,
// columns, bands, ...). The following example demonstrates the setup of an unchecked subvector:

   \code
   blaze::DynamicVector<int> v( 100UL );
   auto sv = subvector( v, 10UL, 20UL, unchecked );  // Creating an unchecked subvector
   \endcode
*/
constexpr Unchecked unchecked;
//*************************************************************************************************

} // namespace blaze

#endif
