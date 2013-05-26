//=================================================================================================
/*!
//  \file blaze/math/CompressedMatrix.h
//  \brief Header file for the complete CompressedMatrix implementation
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

#ifndef _BLAZE_MATH_COMPRESSEDMATRIX_H_
#define _BLAZE_MATH_COMPRESSEDMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/sparse/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/SparseMatrix.h>
#include <blaze/system/Precision.h>
#include <blaze/util/Random.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for CompressedMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of CompressedMatrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
class Rand< CompressedMatrix<Type,SO> >
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand( size_t m, size_t n );
   explicit inline Rand( size_t m, size_t n, size_t nonzeros );
   explicit inline Rand( size_t m, size_t n, Type min, Type max );
   explicit inline Rand( size_t m, size_t n, size_t nonzeros, Type min, Type max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator CompressedMatrix<Type,SO>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CompressedMatrix<Type,SO> matrix_;  //!< The random matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Two argument constructor of the Rand specialization for CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< CompressedMatrix<Type,SO> >::Rand( size_t m, size_t n )
   : matrix_( m, n )  // The random matrix
{
   if( m == 0UL || n == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

   matrix_.reserve( nonzeros );

   while( matrix_.nonZeros() < nonzeros ) {
      matrix_( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Three argument constructor of the Rand specialization for CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \param nonzeros The number of non-zero elements of the random matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< CompressedMatrix<Type,SO> >::Rand( size_t m, size_t n, size_t nonzeros )
   : matrix_( m, n, nonzeros )  // The random matrix
{
   if( nonzeros > m*n )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( m == 0UL || n == 0UL ) return;

   while( matrix_.nonZeros() < nonzeros ) {
      matrix_( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand specialization for CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< CompressedMatrix<Type,SO> >::Rand( size_t m, size_t n, Type min, Type max )
   : matrix_( m, n )  // The random matrix
{
   if( m == 0UL || n == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

   matrix_.reserve( nonzeros );

   while( matrix_.nonZeros() < nonzeros ) {
      matrix_( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Three argument constructor of the Rand specialization for CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \param nonzeros The number of non-zero elements of the random matrix.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< CompressedMatrix<Type,SO> >::Rand( size_t m, size_t n, size_t nonzeros, Type min, Type max )
   : matrix_( m, n, nonzeros )  // The random matrix
{
   if( nonzeros > m*n )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( m == 0UL || n == 0UL ) return;

   while( matrix_.nonZeros() < nonzeros ) {
      matrix_( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random CompressedMatrix.
//
// \return The random matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< CompressedMatrix<Type,SO> >::operator CompressedMatrix<Type,SO>() const
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief MxN single precision matrix.
// \ingroup compressed_matrix
*/
typedef CompressedMatrix<float,false>  CMatMxNf;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief MxN double precision matrix.
// \ingroup compressed_matrix
*/
typedef CompressedMatrix<double,false>  CMatMxNd;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief MxN matrix with system-specific precision.
// \ingroup compressed_matrix
*/
typedef CompressedMatrix<real,false>  CMatMxN;
//*************************************************************************************************

} // namespace blaze

#endif
