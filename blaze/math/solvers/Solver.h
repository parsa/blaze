//=================================================================================================
/*!
//  \file blaze/math/solvers/Solver.h
//  \brief Header file for the base class of all solvers
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

#ifndef _BLAZE_MATH_SOLVERS_SOLVER_H_
#define _BLAZE_MATH_SOLVERS_SOLVER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <limits>
#include <blaze/system/Solvers.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all solver classes.
// \ingroup solver
//
// TODO: description of the Solver class
// TODO: description of its functionality
*/
class Solver
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Solver();
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline size_t getMaxIterations()  const;
   inline size_t getLastIterations() const;
   inline real   getLastPrecision()  const;
   inline real   getThreshold()      const;
   //@}
   //**********************************************************************************************

   //**Set functions***************************************************************************
   /*!\name Set functions */
   //@{
   inline void   setMaxIterations( size_t maxIterations );
   inline void   setThreshold    ( real threshold );
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t maxIterations_;   //!< The maximum number of iterations.
                            /*!< This is the maximum number of iterations the solver will spend
                                 for solving the given problem. */
   size_t lastIterations_;  //!< The number of iterations spent in the last solution process.
   real   lastPrecision_;   //!< The precision of the solution after the solution process.
   real   threshold_;       //!< Precision threshold for the solution.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor.
*/
inline Solver::Solver()
   : maxIterations_ ( solvers::maxIterations )            // The maximum number of iterations
   , lastIterations_( 0 )                                 // The number of iterations spent in the last solution process
   , lastPrecision_ ( std::numeric_limits<real>::max() )  // The precision of the solution after the solution process
   , threshold_     ( solvers::threshold )                // Precision threshold for the solution
{}
//*************************************************************************************************




//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the maximum number of iterations the solver may spend solving the problem.
//
// \return The maximum number of iterations spent in the solver.
*/
inline size_t Solver::getMaxIterations() const
{
   return maxIterations_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of iterations spent in the last solution process.
//
// \return The number of iterations spent in the last solution process.
*/
inline size_t Solver::getLastIterations() const
{
   return lastIterations_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the precision of the solution after the solution process.
//
// \return The precision of the solution after the solution process.
//
// The solver is not enforced to compute the precision after the solution. Instead it can just
// report infinity as the last precision.
*/
inline real Solver::getLastPrecision() const
{
   return lastPrecision_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the threshold that classifies a solution as good enough.
//
// \return The threshold for the solution quality.
*/
inline real Solver::getThreshold() const
{
   return threshold_;
}
//*************************************************************************************************




//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Sets the maximum number of iterations the solver may spend solving the problem.
//
// \param maxIterations The maximum number of iterations spent in the solver.
*/
inline void Solver::setMaxIterations( size_t maxIterations )
{
   maxIterations_ = maxIterations;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets the threshold which classifies a solution as good enough.
//
// \param threshold The threshold for the solution quality.
*/
inline void Solver::setThreshold( real threshold )
{
   threshold_ = threshold;
}
//*************************************************************************************************

} // namespace blaze

#endif
