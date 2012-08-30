//=================================================================================================
/*!
//  \file src/eigen/Daxpy.cpp
//  \brief Source file for the Eigen daxpy product kernel
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
#include <blazemark/eigen/Daxpy.h>
#include <blazemark/eigen/init/Matrix.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Eigen daxpy product kernel.
//
// \param N The size of the vectors for the daxpy product.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the daxpy product by means of the Eigen functionality.
*/
double daxpy( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric_cast;
   using ::Eigen::Dynamic;

   ::blaze::setSeed( seed );

   ::Eigen::Matrix<element_t,Dynamic,1> a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   init( a );

   for( size_t i=0UL; i<N; ++i ) {
      b[i] = element_t(0);
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         b.noalias() += a * element_t(0.001);
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
      std::cerr << " Eigen kernel 'daxpy': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark
