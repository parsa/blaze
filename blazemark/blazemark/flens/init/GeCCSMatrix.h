//=================================================================================================
/*!
//  \file blazemark/flens/init/GeCCSMatrix.h
//  \brief Header file for the FLENS CCS matrix initialization functions
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

#ifndef _BLAZEMARK_FLENS_INIT_GECCSMATRIX_H_
#define _BLAZEMARK_FLENS_INIT_GECCSMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <flens/matrixtypes/general/impl/geccsmatrix.h>
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
void init( ::flens::GeCCSMatrix< ::flens::CCS<Type,::flens::IndexBaseZero<IndexType> > >& m
         , size_t rows, size_t columns, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given CCS matrix.
//
// \param m The CCS matrix to be initialized.
// \param rows The number of rows of the CCS matrix.
// \param columns The number of columns of the CCS matrix.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given CCS matrix with random values. Each row will be filled
// with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type         // Data type of the matrix
        , typename IndexType >  // Index type of the matrix
void init( ::flens::GeCCSMatrix< ::flens::CCS<Type,::flens::IndexBaseZero<IndexType> > >& m
         , size_t rows, size_t columns, size_t nonzeros )
{
   typedef ::flens::IndexBaseZero<IndexType>                              IndexBase;
   typedef ::flens::CoordStorage<Type,::flens::CoordColRowCmp,IndexBase>  Coord;

   ::flens::GeCoordMatrix<Coord> tmp( rows, columns );

   for( IndexType j=tmp.firstCol(); j<=tmp.lastCol(); ++j ) {
      ::blazemark::Indices indices( columns, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         tmp( *it, j ) += ::blaze::rand<Type>( 0, 10 );
      }
   }

   m = tmp;
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark

#endif
