//=================================================================================================
/*!
//  \file src/boost/TDMatSVecMult.cpp
//  \brief Source file for the Boost transpose dense matrix/sparse vector multiplication kernel
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
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <blaze/util/Timing.h>
#include <blazemark/boost/init/CompressedVector.h>
#include <blazemark/boost/init/Matrix.h>
#include <blazemark/boost/TDMatSVecMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace boost {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Boost uBLAS transpose dense matrix/sparse vector multiplication kernel.
//
// \param N The number of rows and columns of the matrix and the size of the vector.
// \param F The number of non-zero elements for the sparse vector.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the transpose dense matrix/sparse vector multiplication by
// means of the Boost uBLAS functionality.
*/
double tdmatsvecmult( size_t N, size_t F, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric::ublas::column_major;

   ::blaze::setSeed( seed );

   ::boost::numeric::ublas::matrix<element_t,column_major> A( N, N );
   ::boost::numeric::ublas::compressed_vector<element_t> a( N );
   ::boost::numeric::ublas::vector<element_t> b( N );
   ::blaze::timing::WcTimer timer;

   init( A );
   init( a, F );

   noalias( b ) = prod( A, a );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         noalias( b ) = prod( A, a );
      }
      timer.end();

      if( b.size() != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Boost uBLAS kernel 'tdmatsvecmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace boost

} // namespace blazemark
