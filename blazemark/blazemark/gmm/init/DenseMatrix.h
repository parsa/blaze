//=================================================================================================
/*!
//  \file blazemark/gmm/init/DenseMatrix.h
//  \brief Header file for the GMM++ dense matrix initialization functions
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

#ifndef _BLAZEMARK_GMM_INIT_DENSEMATRIX_H_
#define _BLAZEMARK_GMM_INIT_DENSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <gmm/gmm_matrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace gmm {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name GMM++ initialization functions */
//@{
template< typename Type >
void init( ::gmm::dense_matrix<Type>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major dense matrix.
//
// \param m The column-major dense matrix to be initialized.
// \return void
//
// This function initializes the given column-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::gmm::dense_matrix<Type>& m )
{
   const size_t M( mat_nrows( m ) );
   const size_t N( mat_ncols( m ) );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark

#endif
