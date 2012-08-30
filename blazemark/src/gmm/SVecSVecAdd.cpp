//=================================================================================================
/*!
//  \file src/gmm/SVecSVecAdd.cpp
//  \brief Source file for the GMM++ sparse vector/sparse vector addition kernel
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
#include <gmm/gmm.h>
#include <blaze/util/Timing.h>
#include <blazemark/gmm/init/RSVector.h>
#include <blazemark/gmm/SVecSVecAdd.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace gmm {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief GMM++ sparse vector/sparse vector addition kernel.
//
// \param N The size of the vectors for the addition.
// \param F The number of non-zero elements for the sparse vectors.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the sparse vector/sparse vector addition by means of the
// GMM++ functionality.
*/
double svecsvecadd( size_t N, size_t F, size_t steps )
{
   using ::blazemark::element_t;

   ::blaze::setSeed( seed );

   ::gmm::rsvector<element_t> a( N ), b( N ), c( N );
   ::blaze::timing::WcTimer timer;

   init( a, F );
   init( b, F );

   ::gmm::add( a, b, c );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         ::gmm::add( a, b, c );
      }
      timer.end();

      if( ::gmm::vect_size(c) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " GMM++ kernel 'svecsvecadd': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark
