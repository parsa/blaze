//=================================================================================================
/*!
//  \file blazemark/eigen/init/SparseMatrix.h
//  \brief Header file for the Eigen sparse matrix initialization functions
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

#ifndef _BLAZEMARK_EIGEN_INIT_SPARSEMATRIX_H_
#define _BLAZEMARK_EIGEN_INIT_SPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <Eigen/Sparse>
#include <blaze/util/Random.h>
#include <blazemark/config/Config.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Indices.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Eigen initialization functions */
//@{
template< typename Type >
void init( ::Eigen::SparseMatrix<Type,::Eigen::RowMajor,EigenSparseIndexType>& m, size_t nonzeros );

template< typename Type >
void init( ::Eigen::SparseMatrix<Type,::Eigen::ColMajor,EigenSparseIndexType>& m, size_t nonzeros );
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
void init( ::Eigen::SparseMatrix<Type,::Eigen::RowMajor,EigenSparseIndexType>& m, size_t nonzeros )
{
   const int M( m.rows() );
   const int N( m.cols() );

   m.reserve( M*nonzeros );

   for( int i=0UL; i<M; ++i ) {
      m.startVec( i );
      ::blazemark::Indices indices( N, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         m.insertBack(i,*it) = ::blaze::rand<Type>();
      }
   }

   m.finalize();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major sparse matrix.
//
// \param m The column-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major sparse matrix with random values. Each column
// will be filled with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::SparseMatrix<Type,::Eigen::ColMajor,EigenSparseIndexType>& m, size_t nonzeros )
{
   const int M( m.rows() );
   const int N( m.cols() );

   m.reserve( N*nonzeros );

   for( int j=0UL; j<N; ++j ) {
      m.startVec( j );
      ::blazemark::Indices indices( M, nonzeros );
      for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         m.insertBack(*it,j) = ::blaze::rand<Type>();
      }
   }

   m.finalize();
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark

#endif
