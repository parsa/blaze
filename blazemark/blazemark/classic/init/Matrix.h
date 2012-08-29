//=================================================================================================
/*!
//  \file blazemark/classic/init/Matrix.h
//  \brief Header file for the classic dense matrix initialization functions
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

#ifndef _BLAZEMARK_CLASSIC_INIT_MATRIX_H_
#define _BLAZEMARK_CLASSIC_INIT_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Random.h>
#include <blazemark/classic/Matrix.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace classic {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Classic initialization functions */
//@{
template< typename Type >
void init( Matrix<Type,false>& m );

template< typename Type >
void init( Matrix<Type,true>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major dense matrix.
//
// \param m The row-major dense matrix to be initialized.
// \return void
//
// This function initializes the given row-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( Matrix<Type,false>& m )
{
   const size_t M( m.rows()    );
   const size_t N( m.columns() );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major dense matrix.
//
// \param m The column-major dense matrix to be initialized.
// \return void
//
// This function initializes the given column-major dense matrix with random values.
*/
template< typename Type >
void init( Matrix<Type,true>& m )
{
   const size_t M( m.rows()    );
   const size_t N( m.columns() );

   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=0UL; i<M; ++i ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************

} // namespace classic

} // namespace blazemark

#endif
