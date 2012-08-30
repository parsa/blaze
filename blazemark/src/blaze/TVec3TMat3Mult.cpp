//=================================================================================================
/*!
//  \file src/blaze/TVec3TMat3Mult.cpp
//  \brief Source file for the Blaze 3D transpose vector/transpose matrix multiplication kernel
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iostream>
#include <vector>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/util/Timing.h>
#include <blazemark/blaze/init/StaticMatrix.h>
#include <blazemark/blaze/init/StaticVector.h>
#include <blazemark/blaze/TVec3TMat3Mult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blaze 3-dimensional transpose vector/transpose matrix multiplication kernel.
//
// \param N The number of 3D vectors to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 3-dimensional transpose vector/transpose matrix
// multiplication by means of the Blaze functionality.
*/
double tvec3tmat3mult( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::rowVector;
   using ::blaze::columnMajor;

   ::blaze::setSeed( seed );

   ::std::vector< ::blaze::StaticVector<element_t,3UL,rowVector> > a( N ), b( N );
   ::std::vector< ::blaze::StaticMatrix<element_t,3UL,3UL,columnMajor> > A( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      init( a[i] );
      init( A[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      b[i] = a[i] * A[i];
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         b[i] = a[i] * A[i];
      }
      timer.end();

      for( size_t i=0UL; i<N; ++i )
         if( b[i][0] < element_t(0) )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blaze kernel 'tvec3tmat3mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark
