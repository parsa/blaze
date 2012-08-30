//=================================================================================================
/*!
//  \file src/blaze/TDVecDVecMult.cpp
//  \brief Source file for the Blaze dense vector/dense vector inner product kernel
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
#include <blaze/math/DynamicVector.h>
#include <blaze/util/Timing.h>
#include <blazemark/blaze/init/DynamicVector.h>
#include <blazemark/blaze/TDVecDVecMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blaze dense vector/dense vector inner product kernel.
//
// \param N The size of the vectors for the inner product.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the dense vector/dense vector inner product by means of
// the Blaze functionality.
*/
double tdvecdvecmult( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::rowVector;
   using ::blaze::columnVector;

   ::blaze::setSeed( seed );

   ::blaze::DynamicVector<element_t,rowVector> a( N );
   ::blaze::DynamicVector<element_t,columnVector> b( N );
   element_t scalar( 0 );
   ::blaze::timing::WcTimer timer;

   init( a );
   init( b );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         scalar += a * b;
      }
      timer.end();

      if( scalar < element_t(0) )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blaze kernel 'tdvecdvecmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark
