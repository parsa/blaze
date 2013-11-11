//=================================================================================================
/*!
//  \file blaze/config/Solvers.h
//  \brief Configuration file for the mathematical solvers
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


namespace blaze {

namespace solvers {

//*************************************************************************************************
/*!\brief Maximum number of iterations of the mathematical solvers.
// \ingroup config
//
// This value specifies the default maximum number of iteration steps the mathematical solvers
// are performing for a single LSE or LCP problem. The choice of the number of iteration steps
// is a consideration between the accuracy of the solution and the time to solution: for a
// large number of iteration steps, the solution will be more accurate, whereas a small number
// will be considerably faster to calculate. For instance, an important property of the PGS
// solver is that the convergence speed of the approximated solution is strongly decreasing
// with a larger number of time steps. In order to noticeably improve the solution, the
// maximum number of time steps has to be considerably increased.
//
// Possible settings for the \a maxIterations variable: \f$ [1..\infty) \f$
*/
const size_t maxIterations = 25000;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Residuum threshold for the LSE/LCP solution of the mathematical solvers.
// \ingroup config
//
// This value specifies the threshold for the residuum calculation during a solution of a
// single LSE or LCP problem. The calculation stops as soon as the approximated solution
// meets the specified threshold value. The choice of the threshold value is a consideration
// between the accuracy of the solution and the time to solution: for a small threshold value
// (e.g. \f$ 1e-9 \f$), the solution will be more accurate, whereas a large threshold value
// (e.g. \f$ 1e-5 \f$) will be considerable faster. For instance, an important property of
// the PGS solver is that the convergance speed of the approximated solution is strongly
// decreasing with a larger number of time steps. The time to solution for a LCP calculation
// using a smaller threshold value will be noticeably larger.
*/
const real threshold = 5E-7;
//*************************************************************************************************

} // namespace solvers

} // namespace pe
