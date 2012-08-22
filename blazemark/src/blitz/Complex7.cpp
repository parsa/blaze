//=================================================================================================
/*!
//  \file src/blitz/Complex7.cpp
//  \brief Source file for the Blitz++ kernel for the complex expression E = ( A + B ) * ( C - D )
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
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/blitz/Complex7.h>
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
/*!\brief Blitz++ kernel for the complex expression E = ( A + B ) * ( C - D ).
//
// \param N The number of rows and columns of the matrices.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the complex expression E = ( A + B ) - ( C - D ) by means of
// the Blitz++ functionality.
*/
double complex7( size_t N, size_t steps )
{
   using ::blazemark::real;
   using ::boost::numeric_cast;

   ::blaze::setSeed( seed );

   ::blitz::Array<real,2> A( N, N, ::blitz::fortranArray );
   ::blitz::Array<real,2> B( N, N, ::blitz::fortranArray );
   ::blitz::Array<real,2> C( N, N, ::blitz::fortranArray );
   ::blitz::Array<real,2> D( N, N, ::blitz::fortranArray );
   ::blitz::Array<real,2> E( N, N, ::blitz::fortranArray );
   ::blitz::firstIndex i;
   ::blitz::secondIndex j;
   ::blitz::thirdIndex k;
   ::blaze::timing::WcTimer timer;

   for( int n=1; n<=static_cast<int>( N ); ++n ) {
      for( int m=1; m<=static_cast<int>( N ); ++m ) {
         A(m,n) = ::blaze::rand<real>();
         B(m,n) = ::blaze::rand<real>();
         C(m,n) = ::blaze::rand<real>();
         D(m,n) = ::blaze::rand<real>();
      }
   }

   {
      ::blitz::Array<real,2> T1( A + B );
      ::blitz::Array<real,2> T2( C - D );
      E = sum( T1(i,k) * T2(k,j), k );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         ::blitz::Array<real,2> T1( A + B );
         ::blitz::Array<real,2> T2( C - D );
         E = sum( T1(i,k) * T2(k,j), k );
      }
      timer.end();

      if( numeric_cast<size_t>( E.rows() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blitz++ kernel 'complex7': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark
