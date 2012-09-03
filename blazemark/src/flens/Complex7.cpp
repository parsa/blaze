//=================================================================================================
/*!
//  \file src/flens/Complex7.cpp
//  \brief Source file for the FLENS kernel for the complex expression E = ( A + B ) * ( C - D )
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
#include <flens/flens.cxx>
#include <blaze/util/Timing.h>
#include <blazemark/flens/Complex7.h>
#include <blazemark/flens/init/GeMatrix.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief FLENS kernel for the complex expression E = ( A + B ) * ( C - D ).
//
// \param N The number of rows and columns of the matrices.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the complex expression E = ( A + B ) * ( C - D ) by means of
// the FLENS functionality.
*/
double complex7( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric_cast;

   typedef ::flens::GeMatrix< ::flens::FullStorage<element_t,::flens::ColMajor> >  ColGeMatrix;

   ::blaze::setSeed( seed );

   ColGeMatrix A( N, N ), B( N, N ), C( N, N ), D( N, N ), E( N, N );
   ::blaze::timing::WcTimer timer;

   init( A );
   init( B );
   init( C );
   init( D );

   {
      const GeMatrix T1( A + B );
      const GeMatrix T2( C - D );
      E = T1 * T2;
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         const GeMatrix T1( A + B );
         const GeMatrix T2( C - D );
         E = T1 * T2;
      }
      timer.end();

      if( numeric_cast<size_t>( E.numRows() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " FLENS kernel 'complex7': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark
