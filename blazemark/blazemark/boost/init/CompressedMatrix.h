//=================================================================================================
/*!
//  \file blazemark/boost/init/CompressedMatrix.h
//  \brief Header file for the Boost sparse matrix initialization functions
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

#ifndef _BLAZEMARK_BOOST_INIT_COMPRESSEDMATRIX_H_
#define _BLAZEMARK_BOOST_INIT_COMPRESSEDMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Indices.h>


namespace blazemark {

namespace boost {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Boost uBLAS initialization functions */
//@{
template< typename Type >
void init( ::boost::numeric::ublas::compressed_matrix<Type,::boost::numeric::ublas::row_major>& m
         , size_t nonzeros );

template< typename Type >
void init( ::boost::numeric::ublas::compressed_matrix<Type,::boost::numeric::ublas::column_major>& m
         , size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major compressed matrix.
//
// \param m The row-major compressed matrix to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given row-major compressed matrix with random values.
// Each row will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::boost::numeric::ublas::compressed_matrix<Type,::boost::numeric::ublas::row_major>& m
         , size_t nonzeros )
{
   const size_t M( m.size1() );
   const size_t N( m.size2() );

   for( size_t i=0UL; i<M; ++i ) {
      ::blazemark::Indices indices( N, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         m(i,*it) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major compressed matrix.
//
// \param m The column-major compressed matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major compressed matrix with random values.
// Each column will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::boost::numeric::ublas::compressed_matrix<Type,::boost::numeric::ublas::column_major>& m
         , size_t nonzeros )
{
   const size_t M( m.size1() );
   const size_t N( m.size2() );

   for( size_t j=0UL; j<N; ++j ) {
      ::blazemark::Indices indices( M, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         m(*it,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************

} // namespace boost

} // namespace blazemark

#endif
