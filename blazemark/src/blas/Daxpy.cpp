//=================================================================================================
/*!
//  \file src/blas/Daxpy.cpp
//  \brief Source file for the BLAS Daxpy product kernel
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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
