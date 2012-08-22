//=================================================================================================
/*!
//  \file src/blitz/TMat3Mat3Mult.cpp
//  \brief Source file for the Blitz++ 3D transpose matrix/matrix multiplication kernel
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
#include <blitz/array.h>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/blitz/TMat3Mat3Mult.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Precision.h>


namespace blazemark {

namespace blitz {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blitz++ 3-dimensional transpose matrix/matrix multiplication kernel.
//
// \param N The number of 3x3 matrices to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 3-dimensional transpose matrix/matrix multiplication by
// means of the Blitz++ functionality.
*/
double tmat3mat3mult( size_t N, size_t steps )
{
   using ::blazemark::real;

   ::blaze::setSeed( seed );

   ::std::vector< ::blitz::Array<real,2> > A( N, ::blitz::Array<real,2>( ::blitz::fortranArray ) );
   ::std::vector< ::blitz::Array<real,2> > B( N );
   ::std::vector< ::blitz::Array<real,2> > C( N, ::blitz::Array<real,2>( ::blitz::fortranArray ) );
   ::blaze::timing::WcTimer timer;

   for( size_t l=0; l<N; ++l ) {
      A[l].resize( 3, 3 );
      C[l].resize( 3, 3 );
      for( int n=1; n<=3; ++n ) {
         for( int m=1; m<=3; ++m ) {
            A[l](m,n) = ::blaze::rand<real>();
         }
      }
   }

   for( size_t l=0; l<N; ++l ) {
      B[l].resize( 3, 3 );
      for( int m=0; m<3; ++m ) {
         for( int n=0; n<3; ++n ) {
            B[l](m,n) = ::blaze::rand<real>();
         }
      }
   }

   for( size_t l=0UL; l<N; ++l ) {
      C[l] = A[l] * B[l];
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, l=0UL; step<steps; ++step, ++l ) {
         if( l == N ) l = 0UL;
         C[l] = A[l] * B[l];
      }
      timer.end();

      for( size_t l=0UL; l<N; ++l )
         if( C[l](0,0) < real(0) )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blitz++ kernel 'tmat3mat3mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark
