//=================================================================================================
/*!
//  \file src/blitz/Complex2.cpp
//  \brief Source file for the Blitz++ kernel for the complex expression d = A * ( a + b + c )
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
#include <blazemark/blitz/Complex2.h>
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
/*!\brief Blitz++ kernel for the complex expression d = A * ( a + b + c ).
//
// \param N The number of rows and columns of the matrix and the size of the vectors.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the complex expression d = A * ( a + b + c ) by means of the
// Blitz++ functionality.
*/
double complex2( size_t N, size_t steps )
{
   using ::blazemark::real;
   using ::boost::numeric_cast;

   ::blaze::setSeed( seed );

   ::blitz::Array<real,2> A( N, N, ::blitz::fortranArray );
   ::blitz::Array<real,1> a( N ), b( N ), c( N ), d( N );
   ::blitz::firstIndex i;
   ::blitz::secondIndex j;
   ::blaze::timing::WcTimer timer;

   for( int n=1; n<=static_cast<int>( N ); ++n ) {
      for( int m=1; m<=static_cast<int>( N ); ++m ) {
         A(m,n) = ::blaze::rand<real>();
      }
   }

   for( int m=0; m<static_cast<int>( N ); ++m ) {
      a(m) = ::blaze::rand<real>();
      b(m) = ::blaze::rand<real>();
      c(m) = ::blaze::rand<real>();
   }

   {
      ::blitz::Array<real,1> tmp( a + b + c );
      d = sum( A(i,j) * tmp(j), j );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         ::blitz::Array<real,1> tmp( a + b + c );
         d = sum( A(i,j) * tmp(j), j );
      }
      timer.end();

      if( numeric_cast<size_t>( d.size() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blitz++ kernel 'complex2': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark
