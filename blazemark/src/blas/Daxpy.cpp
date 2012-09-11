//=================================================================================================
/*!
//  \file src/blas/Daxpy.cpp
//  \brief Source file for the BLAS Daxpy product kernel
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
#include <blaze/math/DynamicVector.h>
#include <blaze/util/Timing.h>
#include <blazemark/blas/init/DynamicVector.h>
#include <blazemark/blas/Daxpy.h>
#include <blazemark/system/BLAS.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blas {

//=================================================================================================
//
//  COMPUTE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Kernel function for single precision vectors.
//
// \param N Size of the vectors.
// \param alpha Scalar factor for \f$ \alpha op(x) \f$.
// \param X Pointer to the first element of vector X.
// \param incX Use every incX'th element of vector X.
// \param Y Pointer to the first element of vector Y.
// \param incY Use every incY'th element of vector Y.
// \return void
*/
inline void daxpy( const int N, const float alpha,
                   const float *X, const int incX, float *Y, const int incY )
{
   cblas_saxpy( N, alpha, X, incX, Y, incY );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Kernel function for double precision vectors.
//
// \param N Size of the vectors.
// \param alpha Scalar factor for \f$ \alpha op(x) \f$.
// \param X Pointer to the first element of vector X.
// \param incX Use every incX'th element of vector X.
// \param Y Pointer to the first element of vector Y.
// \param incY Use every incY'th element of vector Y.
// \return void
*/
inline void daxpy( const int N, const double alpha,
                   const double *X, const int incX, double *Y, const int incY )
{
   cblas_daxpy( N, alpha, X, incX, Y, incY );
}
//*************************************************************************************************




//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief BLAS daxpy product kernel.
//
// \param N The size of the vectors for the daxpy product.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the daxpy product by means of BLAS functionality.
*/
double daxpy( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::columnVector;

   ::blaze::setSeed( seed );

   ::blaze::DynamicVector<element_t,columnVector> a( N ), b( N, 0 );
   ::blaze::timing::WcTimer timer;

   init( a );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         daxpy( N, element_t(3), a.data(), 1, b.data(), 1 );
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
      std::cerr << " Blaze kernel 'daxpy': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blas

} // namespace blazemark
