//=================================================================================================
/*!
//  \file blaze/math/DecompositionFlag.h
//  \brief Header file for the dense matrix decomposition flags
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

#ifndef _BLAZE_MATH_DECOMPOSITIONFLAG_H_
#define _BLAZE_MATH_DECOMPOSITIONFLAG_H_


namespace blaze {

//=================================================================================================
//
//  DECOMPOSITION FLAG VALUES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Decomposition flag.
// \ingroup dense_matrix
//
// The DecompositionFlag type enumeration represents the different types of matrix decomposition
// that are available within the Blaze library. The following flags are available:
//
//  - \a byLU: The default decomposition algorithm for general square matrices. It decomposes a
//             matrix into a lower unitriangular matrix \c L, an upper triangular matrix \c U,
//             and a permutation matrix \c P (\f$ A = L U P \f$). If no permutations are required,
//             \c P is the identity matrix.
//  - \a byCholesky: An optimized decomposition for symmetric positive definite matrices. It
//             decomposes a given matrix into either \f$ A = L^T L \f$, where \c L is a lower
//             triangular matrix, or \f$ A = U^T U \f$, where \c U is an upper triangular matrix.
//  - \a byQR: A very general decomposition for M-by-N matrices. It decomposes the matrix into
//             \f$ A = Q R \f$, where \c Q is an orthogonal matrix of size M-by-M, and \c R is
//             an upper triangular matrix of size M-by-N.
*/
enum DecompositionFlag
{
   byLU       = 0,  //!< Flag for the LU decomposition.
   byCholesky = 1,  //!< Flag for the Cholesky decomposition.
   byQR       = 2   //!< Flag for the QR decomposition.
};
//*************************************************************************************************

} // namespace blaze

#endif
