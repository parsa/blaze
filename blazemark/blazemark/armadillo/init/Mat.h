//=================================================================================================
/*!
//  \file blazemark/armadillo/init/Mat.h
//  \brief Header file for the Armadillo dense matrix initialization functions
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

#ifndef _BLAZEMARK_ARMADILLO_INIT_MAT_H_
#define _BLAZEMARK_ARMADILLO_INIT_MAT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <armadillo>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace armadillo {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Armadillo initialization functions */
//@{
template< typename Type >
void init( ::arma::Mat<Type>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major dense matrix.
//
// \param m The dense column-major dense matrix to be initialized.
// \return void
//
// This function initializes the given dense column-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::arma::Mat<Type>& m )
{
   const size_t M( m.n_rows );
   const size_t N( m.n_cols );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************

} // namespace armadillo

} // namespace blazemark

#endif
