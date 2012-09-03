//=================================================================================================
/*!
//  \file blazemark/flens/init/GeMatrix.h
//  \brief Header file for the FLENS dense matrix initialization functions
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

#ifndef _BLAZEMARK_FLENS_INIT_GEMATRIX_H_
#define _BLAZEMARK_FLENS_INIT_GEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name FLENS initialization functions */
//@{
template< typename Type >
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::RowMajor> >& m );

template< typename Type >
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::ColMajor> >& m );
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
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::RowMajor> >& m )
{
   typedef typename ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::RowMajor> >::IndexType  IndexType;

   for( IndexType i=m.firstRow(); i<=m.lastRow(); ++i ) {
      for( IndexType j=m.firstCol(); j<=m.lastCol(); ++j ) {
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
template< typename Type >  // Data type of the matrix
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::ColMajor> >& m )
{
   typedef typename ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::ColMajor> >::IndexType  IndexType;

   for( IndexType j=m.firstCol(); j<=m.lastCol(); ++j ) {
      for( IndexType i=m.firstRow(); i<=m.lastRow(); ++i ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark

#endif
