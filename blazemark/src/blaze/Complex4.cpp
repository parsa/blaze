//=================================================================================================
/*!
//  \file src/blaze/Complex4.cpp
//  \brief Source file for the Blaze kernel for the complex expression b += s * A * a
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/Timing.h>
#include <blazemark/blaze/Complex4.h>
#include <blazemark/blaze/init/DynamicMatrix.h>
#include <blazemark/blaze/init/DynamicVector.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blaze kernel for the complex expression b += s * A * a.
//
// \param N The number of rows and columns of the matrices and the size of the vectors.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the complex expression b += s * A * a by means of the
// Blaze functionality.
*/
double complex4( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::columnVector;
   using ::blaze::columnMajor;

   ::blaze::setSeed( seed );

   ::blaze::DynamicMatrix<element_t,columnMajor> A( N, N );
   ::blaze::DynamicVector<element_t,columnVector> a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   init( A );
   init( a );

   b = element_t(0);

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         b += element_t(2.2) * A * a;
      }
      timer.end();

      if( b.size() != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blaze kernel 'complex4': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark
