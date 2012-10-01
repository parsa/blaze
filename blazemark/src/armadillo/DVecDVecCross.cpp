//=================================================================================================
/*!
//  \file src/armadillo/DVecDVecCross.cpp
//  \brief Source file for the Armadillo dense vector/dense vector cross product kernel
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
#include <armadillo>
#include <blaze/util/Timing.h>
#include <blazemark/armadillo/DVecDVecCross.h>
#include <blazemark/armadillo/init/Col.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace armadillo {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Armadillo dense vector/dense vector cross product kernel.
//
// \param N The number of cross products to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the dense vector/dense vector cross product by means of
// the Armadillo functionality.
*/
double dvecdveccross( size_t N, size_t steps )
{
   using ::blazemark::element_t;

   ::blaze::setSeed( seed );

   ::std::vector< ::arma::Col<element_t> > a( N ), b( N ), c( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0; i<N; ++i ) {
      a[i].resize( 3UL );
      b[i].resize( 3UL );
      c[i].resize( 3UL );
      init( a[i] );
      init( b[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      c[i] = cross( a[i], b[i] );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         c[i] = cross( a[i], b[i] );
      }
      timer.end();

      for( size_t i=0UL; i<N; ++i )
         if( c[i].n_elem != 3UL )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Armadillo kernel 'dvecdveccross': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace armadillo

} // namespace blazemark
