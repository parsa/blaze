//=================================================================================================
/*!
//  \file src/clike/Mat3Vec3Mult.cpp
//  \brief Source file for the C-like 3D matrix/vector multiplication kernel
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iostream>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/clike/init/Matrix.h>
#include <blazemark/clike/init/Vector.h>
#include <blazemark/clike/Matrix.h>
#include <blazemark/clike/Mat3Vec3Mult.h>
#include <blazemark/clike/Vector.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace clike {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief C-like 3-dimensional matrix/vector multiplication kernel.
//
// \param N The number of 3D vectors to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 3-dimensional matrix/vector multiplication by means of a
// C-like implementation.
*/
double mat3vec3mult( size_t N, size_t steps )
{
   using ::blazemark::element_t;

   ::blaze::setSeed( seed );

   ::std::vector< Matrix<element_t,3UL,3UL> > A( N );
   ::std::vector< Vector<element_t,3UL> > a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      init( A[i] );
      init( a[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      b[i].v[0] = A[i].v[0][0] * a[i].v[0] + A[i].v[0][1] * a[i].v[1] + A[i].v[0][2] * a[i].v[2];
      b[i].v[1] = A[i].v[1][0] * a[i].v[0] + A[i].v[1][1] * a[i].v[1] + A[i].v[1][2] * a[i].v[2];
      b[i].v[2] = A[i].v[2][0] * a[i].v[0] + A[i].v[2][1] * a[i].v[1] + A[i].v[2][2] * a[i].v[2];
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         b[i].v[0] = A[i].v[0][0] * a[i].v[0] + A[i].v[0][1] * a[i].v[1] + A[i].v[0][2] * a[i].v[2];
         b[i].v[1] = A[i].v[1][0] * a[i].v[0] + A[i].v[1][1] * a[i].v[1] + A[i].v[1][2] * a[i].v[2];
         b[i].v[2] = A[i].v[2][0] * a[i].v[0] + A[i].v[2][1] * a[i].v[1] + A[i].v[2][2] * a[i].v[2];
      }
      timer.end();

      for( size_t i=0UL; i<N; ++i )
         if( b[i].v[0] < element_t(0) )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " C-like kernel 'mat3vec3mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace clike

} // namespace blazemark
