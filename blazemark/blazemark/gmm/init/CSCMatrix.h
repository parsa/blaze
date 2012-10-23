//=================================================================================================
/*!
//  \file blazemark/gmm/init/CSCMatrix.h
//  \brief Header file for the GMM++ column-major sparse matrix initialization functions
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

#ifndef _BLAZEMARK_GMM_INIT_CSCMATRIX_H_
#define _BLAZEMARK_GMM_INIT_CSCMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <gmm/gmm_matrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Indices.h>


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
void init( ::gmm::csc_matrix<Type>& m, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major sparse matrix.
//
// \param m The column-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major sparse matrix with random values.
// Each column will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::gmm::csc_matrix<Type>& m, size_t nonzeros )
{
   const size_t M( mat_nrows(m) );
   const size_t N( mat_ncols(m) );

   ::gmm::col_matrix< ::gmm::wsvector<Type> > tmp( M, N );

   if( structure == band )
   {
      const size_t drange( nonzeros / 2UL );
      const size_t urange( ( nonzeros % 2UL )?( drange ):( drange-1UL ) );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( j >= urange )?( j-urange ):( 0UL ) );
         const size_t iend  ( ( j+drange+1UL < M )?( j+drange+1UL ):( M ) );

         for( size_t i=ibegin; i<iend; ++i ) {
            tmp(i,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<N; ++j ) {
         ::blazemark::Indices indices( M, nonzeros );
         for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
            tmp(*it,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }

   copy( tmp, m );
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark

#endif
