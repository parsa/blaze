//=================================================================================================
/*!
//  \file src/eigen/CG.cpp
//  \brief Source file for the Eigen conjugate gradient kernel
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
#include <Eigen/Sparse>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/eigen/CG.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Precision.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Eigen conjugate gradient kernel.
//
// \param N The number of rows and columns of the 2D discretized grid.
// \param steps The number of solving steps to perform.
// \param iterations The number of iterations to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the conjugate gradient method by means of the Eigen
// functionality.
*/
double cg( size_t N, size_t steps, size_t iterations )
{
   using ::blazemark::real;
   using ::boost::numeric_cast;
   using ::Eigen::Dynamic;
   using ::Eigen::RowMajor;

   ::blaze::setSeed( seed );

   const size_t NN( N*N );

   ::Eigen::SparseMatrix<real,RowMajor,EigenSparseIndexType> A( NN, NN );
   ::Eigen::Matrix<real,Dynamic,1> x( NN ), b( NN ), r( NN ), d( NN ), h( NN ), init( NN );
   real alpha, beta, delta;
   ::blaze::timing::WcTimer timer;

   A.reserve( NN*5UL );

   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         A.startVec( i*N+j );
         if( i > 0UL   ) A.insertBack(i*N+j,(i-1UL)*N+j) = -1.0;  // Top neighbor
         if( j > 0UL   ) A.insertBack(i*N+j,i*N+j-1UL  ) = -1.0;  // Left neighbor
         A.insertBack(i*N+j,i*N+j) = 4.0;
         if( j < N-1UL ) A.insertBack(i*N+j,i*N+j+1UL  ) = -1.0;  // Right neighbor
         if( i < N-1UL ) A.insertBack(i*N+j,(i+1UL)*N+j) = -1.0;  // Bottom neighbor
      }
   }

   A.finalize();

   for( size_t i=0UL; i<NN; ++i ) {
      b[i]    = real(0);
      init[i] = ::blaze::rand<real>();
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step )
      {
         x.noalias() = init;
         r.noalias() = A * x - b;
         delta = r.transpose() * r;
         d.noalias() = -r;

         for( size_t iteration=0UL; iteration<iterations; ++iteration )
         {
            h = A * d;
            alpha = delta / ( d.transpose() * h );
            x += alpha * d;
            r += alpha * h;
            beta = r.transpose() * r;
            d.noalias() = ( beta / delta ) * d - r;
            delta = beta;
         }
      }

      timer.end();

      if( numeric_cast<size_t>( x.size() ) != NN )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Eigen kernel 'cg': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark
