//=================================================================================================
/*!
//  \file src/boost/Complex3.cpp
//  \brief Source file for the Boost kernel for the complex expression c = A * B * ( a + b )
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
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <blaze/util/Timing.h>
#include <blazemark/boost/Complex3.h>
#include <blazemark/boost/init/Matrix.h>
#include <blazemark/boost/init/Vector.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace boost {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Boost uBLAS kernel for the complex expression c = A * B * ( a + b ).
//
// \param N The number of rows and columns of the matrices and the size of the vectors.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the complex expression c = A * B * ( a + b ) by means of the
// Boost uBLAS functionality.
*/
double complex3( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric::ublas::column_major;

   ::blaze::setSeed( seed );

   ::boost::numeric::ublas::matrix<element_t,column_major> A( N, N ), B( N, N );
   ::boost::numeric::ublas::vector<element_t> a( N ), b( N ), c( N );
   ::blaze::timing::WcTimer timer;

   init( A );
   init( B );
   init( a );
   init( b );

   {
      ::boost::numeric::ublas::vector<element_t> tmp( prod( B, ( a + b ) ) );
      noalias( c )   = prod( A, tmp );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         ::boost::numeric::ublas::vector<element_t> tmp( prod( B, ( a + b ) ) );
         noalias( c )   = prod( A, tmp );
      }
      timer.end();

      if( c.size() != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Boost uBLAS kernel 'complex3': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace boost

} // namespace blazemark
