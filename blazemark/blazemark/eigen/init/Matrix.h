//=================================================================================================
/*!
//  \file blazemark/eigen/init/Matrix.h
//  \brief Header file for the Eigen dense matrix initialization functions
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

#ifndef _BLAZEMARK_EIGEN_INIT_MATRIX_H_
#define _BLAZEMARK_EIGEN_INIT_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <Eigen/Dense>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


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
void init( ::Eigen::Matrix<Type,3,1>& v );

template< typename Type >
void init( ::Eigen::Matrix<Type,6,1>& v );

template< typename Type >
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,1>& v );

template< typename Type >
void init( ::Eigen::Matrix<Type,3,3,::Eigen::RowMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,3,3,::Eigen::ColMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,6,6,::Eigen::RowMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,6,6,::Eigen::ColMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::RowMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::ColMajor>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given 3D vector.
//
// \param v The 3D vector to be initialized.
// \return void
//
// This function initializes the given 3D vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::Eigen::Matrix<Type,3,1>& v )
{
   for( int i=0; i<3; ++i ) {
      v[i] = ::blaze::rand<Type>();
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given 6D vector.
//
// \param v The 6D vector to be initialized.
// \return void
//
// This function initializes the given 6D vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::Eigen::Matrix<Type,6,1>& v )
{
   for( int i=0; i<6; ++i ) {
      v[i] = ::blaze::rand<Type>();
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given dense vector.
//
// \param v The dense vector to be initialized.
// \return void
//
// This function initializes the given dense vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,1>& v )
{
   const int N( v.size() );

   for( int i=0; i<N; ++i ) {
      v[i] = ::blaze::rand<Type>();
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major 3x3 matrix.
//
// \param m The row-major 3x3 matrix to be initialized.
// \return void
//
// This function initializes the given row-major 3x3 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,3,3,::Eigen::RowMajor>& m )
{
   for( int i=0; i<3; ++i ) {
      for( int j=0; j<3; ++j ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major 3x3 matrix.
//
// \param m The column-major 3x3 matrix to be initialized.
// \return void
//
// This function initializes the given column-major 3x3 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,3,3,::Eigen::ColMajor>& m )
{
   for( int i=0; i<3; ++i ) {
      for( int j=0; j<3; ++j ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major 6x6 matrix.
//
// \param m The row-major 6x6 matrix to be initialized.
// \return void
//
// This function initializes the given row-major 6x6 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,6,6,::Eigen::RowMajor>& m )
{
   for( int i=0; i<6; ++i ) {
      for( int j=0; j<6; ++j ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major 6x6 matrix.
//
// \param m The column-major 6x6 matrix to be initialized.
// \return void
//
// This function initializes the given column-major 6x6 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,6,6,::Eigen::ColMajor>& m )
{
   for( int i=0; i<6; ++i ) {
      for( int j=0; j<6; ++j ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
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
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::RowMajor>& m )
{
   const int M( m.rows() );
   const int N( m.cols() );

   for( int i=0; i<M; ++i ) {
      for( int j=0; j<N; ++j ) {
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
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::ColMajor>& m )
{
   const int M( m.rows() );
   const int N( m.cols() );

   for( int j=0; j<N; ++j ) {
      for( int i=0; i<M; ++i ) {
         m(i,j) = ::blaze::rand<Type>();
      }
   }
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark

#endif
