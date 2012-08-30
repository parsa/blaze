//=================================================================================================
/*!
//  \file src/blitz/DMatDMatMult.cpp
//  \brief Source file for the Blitz++ dense matrix/dense matrix multiplication kernel
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
#include <blitz/array.h>
#include <boost/cast.hpp>
#include <blaze/util/Timing.h>
#include <blazemark/blitz/DMatDMatMult.h>
#include <blazemark/blitz/init/Array.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blitz {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blitz++ dense matrix/dense matrix multiplication kernel.
//
// \param N The number of rows and columns of the matrices.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the dense matrix/dense matrix multiplication by means of
// the Blitz++ functionality.
*/
double dmatdmatmult( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric_cast;

   ::blaze::setSeed( seed );

   ::blitz::Array<element_t,2> A( N, N ), B( N, N ), C( N, N );
   ::blitz::firstIndex i;
   ::blitz::secondIndex j;
   ::blitz::thirdIndex k;
   ::blaze::timing::WcTimer timer;

   initRowMajorMatrix( A );
   initRowMajorMatrix( B );

   C = sum( A(i,k) * B(k,j), k );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         C = sum( A(i,k) * B(k,j), k );
      }
      timer.end();

      if( numeric_cast<size_t>( C.rows() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blitz++ kernel 'dmatdmatmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark
