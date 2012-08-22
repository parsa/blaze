//=================================================================================================
/*!
//  \file src/boost/TMat3Mat3Mult.cpp
//  \brief Source file for the Boost 3D transpose matrix/matrix multiplication kernel
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
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/boost/TMat3Mat3Mult.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Precision.h>


namespace blazemark {

namespace boost {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Boost uBLAS 3-dimensional transpose matrix/matrix multiplication kernel.
//
// \param N The number of 3x3 matrices to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 3-dimensional transpose matrix/matrix multiplication by
// means of the Boost uBLAS functionality.
*/
double tmat3mat3mult( size_t N, size_t steps )
{
   using ::blazemark::real;
   using ::boost::numeric::ublas::row_major;
   using ::boost::numeric::ublas::column_major;

   ::blaze::setSeed( seed );

   ::std::vector< ::boost::numeric::ublas::matrix<real,column_major> > A( N ), C( N );
   ::std::vector< ::boost::numeric::ublas::matrix<real,row_major> > B( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      A[i].resize( 3UL, 3UL );
      B[i].resize( 3UL, 3UL );
      C[i].resize( 3UL, 3UL );
      for( size_t k=0UL; k<3UL; ++k ) {
         for( size_t j=0UL; j<3UL; ++j ) {
            A[i](j,k) = ::blaze::rand<real>();
            B[i](j,k) = ::blaze::rand<real>();
         }
      }
   }

   for( size_t i=0UL; i<N; ++i ) {
      noalias( C[i] ) = prod( A[i], B[i] );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         noalias( C[i] ) = prod( A[i], B[i] );
      }
      timer.end();

      for( size_t i=0UL; i<N; ++i )
         if( C[i](0,0) < real(0) )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Boost uBLAS kernel 'tmat3mat3mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace boost

} // namespace blazemark
