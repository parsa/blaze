//=================================================================================================
/*!
//  \file src/eigen/SMatSMatAdd.cpp
//  \brief Source file for the Eigen sparse matrix/sparse matrix addition kernel
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
#include <Eigen/Sparse>
#include <blaze/util/Timing.h>
#include <blazemark/eigen/init/SparseMatrix.h>
#include <blazemark/eigen/SMatSMatAdd.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Eigen sparse matrix/sparse matrix addition kernel.
//
// \param N The number of rows and columns of the matrices.
// \param F The number of non-zero elements in each row of the sparse matrices.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the sparse matrix/sparse matrix addition by means of the
// Eigen functionality.
*/
double smatsmatadd( size_t N, size_t F, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric_cast;
   using ::Eigen::RowMajor;

   ::blaze::setSeed( seed );

   ::Eigen::SparseMatrix<element_t,RowMajor,EigenSparseIndexType> A( N, N ), B( N, N ), C( N, N );
   ::blaze::timing::WcTimer timer;

   init( A, F );
   init( B, F );

   C = A + B;

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         C = A + B;
      }
      timer.end();

      if( numeric_cast<size_t>( C.rows() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Eigen kernel 'smatsmatadd': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark
