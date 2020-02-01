//=================================================================================================
/*!
//  \file blaze/math/ReductionFlag.h
//  \brief Header file for the reduction flag enumeration
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_REDUCTIONFLAG_H_
#define _BLAZE_MATH_REDUCTIONFLAG_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  REDUCTION FLAGS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Reduction flag for row-wise or column-wise reduction operations.
// \ingroup math
//
// This flag can be used to perform row-wise or column-wise reduction operations on matrices.
// The following example shows both the row-wise and column-wise summation of a dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicMatrix<int,rowMajor> A{ { 4, 1, 2 }, { -2, 0, 3 } };
   blaze::DynamicVector<int,columnVector> v1;
   blaze::DynamicVector<int,rowVector> v2;

   v1 = sum<rowwise>( A );     // Results in ( 7, 1 )
   v2 = sum<columnwise>( A );  // Results in ( 2, 1, 5 )
   \endcode
*/
enum ReductionFlag : size_t
{
   columnwise = 0UL,  //!< Flag for column-wise reduction operations.
   rowwise    = 1UL   //!< Flag for row-wise reduction operations.
};
//*************************************************************************************************

} // namespace blaze

#endif
