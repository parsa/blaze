//=================================================================================================
/*!
//  \file src/clike/SMatDVecMult.cpp
//  \brief Source file for the C-like sparse matrix/dense vector multiplication kernel
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
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/clike/SMatDVecMult.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Precision.h>
#include <blazemark/util/Indices.h>


namespace blazemark {

namespace clike {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief C-like sparse matrix/dense vector multiplication kernel.
//
// \param N The number of rows and columns of the matrix and the size of the vector.
// \param F The number of non-zero elements in each row of the sparse matrix.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the sparse matrix/dense vector multiplication by means of
// a C-like implementation.
*/
double smatdvecmult( size_t N, size_t F, size_t steps )
{
   using ::blazemark::real;

   ::blaze::setSeed( seed );

   real* value = new real[F*N];
   size_t* index = new size_t[F*N];
   size_t* row = new size_t[N+1UL];
   real* a = new real[N];
   real* b = new real[N];
   ::blaze::timing::WcTimer timer;
   size_t counter( 0 );

   for( size_t i=0UL; i<N; ++i ) {
      ::blazemark::Indices indices( N, F );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         value[counter] = ::blaze::rand<real>();
         index[counter] = *it;
         ++counter;
      }
      row[i] = i*F;
   }
   row[N] = N*F;

   for( size_t i=0UL; i<N; ++i ) {
      a[i] = ::blaze::rand<real>();
   }

   for( size_t i=0UL; i<N; ++i ) {
      b[i] = real(0);
      const size_t begin( row[i    ] );
      const size_t end  ( row[i+1UL] );
      for( size_t j=begin; j!=end; ++j ) {
         b[i] += value[j] * a[index[j]];
      }
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         for( size_t i=0UL; i<N; ++i ) {
            b[i] = real(0);
            const size_t begin( row[i    ] );
            const size_t end  ( row[i+1UL] );
            for( size_t j=begin; j!=end; ++j ) {
               b[i] += value[j] * a[index[j]];
            }
         }
      }
      timer.end();

      if( b[0] < real(0) )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   delete[] value;
   delete[] index;
   delete[] row;
   delete[] a;
   delete[] b;

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " C-like kernel 'smatdvecmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace clike

} // namespace blazemark
