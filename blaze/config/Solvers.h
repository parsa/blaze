//=================================================================================================
/*!
//  \file blaze/config/Solvers.h
//  \brief Configuration file for the mathematical solvers
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
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
