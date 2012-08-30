//=================================================================================================
/*!
//  \file src/eigen/TDMatDVecMult.cpp
//  \brief Source file for the Eigen transpose dense matrix/dense vector multiplication kernel
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
#include <boost/cast.hpp>
#include <Eigen/Dense>
#include <blaze/util/Timing.h>
#include <blazemark/eigen/init/Matrix.h>
#include <blazemark/eigen/TDMatDVecMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Eigen transpose dense matrix/dense vector multiplication kernel.
//
// \param N The number of rows and columns of the matrix and the size of the vector.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the transpose dense matrix/dense vector multiplication by
// means of the Eigen functionality.
*/
double tdmatdvecmult( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric_cast;
   using ::Eigen::Dynamic;
   using ::Eigen::ColMajor;

   ::blaze::setSeed( seed );

   ::Eigen::Matrix<element_t,Dynamic,Dynamic,ColMajor> A( N, N );
   ::Eigen::Matrix<element_t,Dynamic,1> a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   init( A );
   init( a );

   b.noalias() = A * a;

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         b.noalias() = A * a;
      }
      timer.end();

      if( numeric_cast<size_t>( b.size() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Eigen kernel 'tdmatdvecmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark
