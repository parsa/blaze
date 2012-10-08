//=================================================================================================
/*!
//  \file blazemark/blaze/init/StaticMatrix.h
//  \brief Header file for the Blaze static matrix initialization functions
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

#ifndef _BLAZEMARK_BLAZE_INIT_STATICMATRIX_H_
#define _BLAZEMARK_BLAZE_INIT_STATICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Blaze initialization functions */
//@{
template< typename Type, size_t M, size_t N >
void init( ::blaze::StaticMatrix<Type,M,N,::blaze::rowMajor>& m );

template< typename Type, size_t M, size_t N >
void init( ::blaze::StaticMatrix<Type,M,N,::blaze::columnMajor>& m );

template< typename Type, size_t M, size_t N, bool SO >
void init( ::std::vector< ::blaze::StaticMatrix<Type,M,N,SO> >& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major static matrix.
//
// \param m The row-major static matrix to be initialized.
// \return void
//
// This function initializes the given row-major static matrix with random values.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
void init( ::blaze::StaticMatrix<Type,M,N,::blaze::rowMajor>& m )
{
   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major static matrix.
//
// \param m The column-major static matrix to be initialized.
// \return void
//
// This function initializes the given column-major static matrix with random values.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N >     // Number of columns
void init( ::blaze::StaticMatrix<Type,M,N,::blaze::columnMajor>& m )
{
   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given vector of static matrices.
//
// \param v The vector of static matrices to be initialized.
// \return void
//
// This function initializes all static matrices in the given vector with random values.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
void init( ::std::vector< ::blaze::StaticMatrix<Type,M,N,SO> >& v )
{
   const size_t size( v.size() );

   for( size_t i=0UL; i<size; ++i ) {
      init( v[i] );
   }
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark

#endif
