//=================================================================================================
/*!
//  \file src/blaze/SMatScalarMult.cpp
//  \brief Source file for the Blaze sparse matrix/scalar multiplication kernel
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/util/Timing.h>
#include <blazemark/blaze/init/CompressedMatrix.h>
#include <blazemark/blaze/SMatScalarMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blaze sparse matrix/scalar multiplication kernel.
//
// \param N The number of rows and columns of the matrix.
// \param F The number of non-zero elements in each row of the sparse matrix.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the sparse matrix/scalar multiplication by means of the
// Blaze functionality.
*/
double smatscalarmult( size_t N, size_t F, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::rowMajor;

   ::blaze::setSeed( seed );

   ::blaze::CompressedMatrix<element_t,rowMajor> A( N, N, N*F ), B( N, N );
   ::blaze::timing::WcTimer timer;

   init( A, F );

   B = A * element_t(2.2);

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         B = A * element_t(2.2);
      }
      timer.end();

      if( B.rows() != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blaze kernel 'smatscalarmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark
