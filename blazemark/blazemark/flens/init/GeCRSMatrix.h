//=================================================================================================
/*!
//  \file blazemark/flens/init/GeCRSMatrix.h
//  \brief Header file for the FLENS CRS matrix initialization functions
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

#ifndef _BLAZEMARK_FLENS_INIT_GECRSMATRIX_H_
#define _BLAZEMARK_FLENS_INIT_GECRSMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <flens/matrixtypes/general/impl/gecrsmatrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Indices.h>


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
template< typename Type, typename IndexType >
void init( ::flens::GeCRSMatrix< ::flens::CRS<Type,::flens::IndexBaseZero<IndexType> > >& m
         , size_t rows, size_t columns, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given CRS matrix.
//
// \param m The CRS matrix to be initialized.
// \param rows The number of rows of the CRS matrix.
// \param columns The number of columns of the CRS matrix.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given CRS matrix with random values. Each row will be filled
// with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type         // Data type of the matrix
        , typename IndexType >  // Index type of the matrix
void init( ::flens::GeCRSMatrix< ::flens::CRS<Type,::flens::IndexBaseZero<IndexType> > >& m
         , size_t rows, size_t columns, size_t nonzeros )
{
   typedef ::flens::IndexBaseZero<IndexType>                              IndexBase;
   typedef ::flens::CoordStorage<Type,::flens::CoordRowColCmp,IndexBase>  Coord;

   ::flens::GeCoordMatrix<Coord> tmp( rows, columns );

   for( IndexType i=tmp.firstRow(); i<=tmp.lastRow(); ++i ) {
      ::blazemark::Indices indices( columns, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         tmp( i, *it ) += ::blaze::rand<Type>( 0, 10 );
      }
   }

   m = tmp;
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark

#endif
