//=================================================================================================
/*!
//  \file blazemark/mtl/init/Compressed2D.h
//  \brief Header file for the MTL compressed matrix initialization functions
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

#ifndef _BLAZEMARK_MTL_INIT_COMPRESSED2D_H_
#define _BLAZEMARK_MTL_INIT_COMPRESSED2D_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Indices.h>


namespace blazemark {

namespace mtl {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name MTL initialization functions */
//@{
template< typename Type >
void init( ::mtl::compressed2D< Type, ::mtl::matrix::parameters< ::mtl::tag::row_major > >& m
         , size_t nonzeros );

template< typename Type >
void init( ::mtl::compressed2D< Type, ::mtl::matrix::parameters< ::mtl::tag::col_major > >& m
         , size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major sparse matrix.
//
// \param m The row-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given row-major sparse matrix with random values. Each row
// will be filled with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::mtl::compressed2D< Type, ::mtl::matrix::parameters< ::mtl::tag::row_major > >& m
         , size_t nonzeros )
{
   typedef ::mtl::tag::row_major                      row_major;
   typedef ::mtl::matrix::parameters<row_major>       row_parameters;
   typedef ::mtl::compressed2D<Type,row_parameters>   row_compressed2D;
   typedef ::mtl::matrix::inserter<row_compressed2D>  row_inserter;

   const size_t M( num_rows( m ) );
   const size_t N( num_cols( m ) );

   row_inserter ins( m );

   for( size_t i=0UL; i<M; ++i ) {
      ::blazemark::Indices indices( N, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         ins[i][*it] = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
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
void init( ::mtl::compressed2D< Type, ::mtl::matrix::parameters< ::mtl::tag::col_major > >& m
         , size_t nonzeros )
{
   typedef ::mtl::tag::col_major                      col_major;
   typedef ::mtl::matrix::parameters<col_major>       col_parameters;
   typedef ::mtl::compressed2D<Type,col_parameters>   col_compressed2D;
   typedef ::mtl::matrix::inserter<col_compressed2D>  col_inserter;

   const size_t M( num_rows( m ) );
   const size_t N( num_cols( m ) );

   col_inserter ins( m );

   for( size_t j=0UL; j<N; ++j ) {
      ::blazemark::Indices indices( M, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         ins[*it][j] = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************

} // namespace mtl

} // namespace blazemark

#endif
