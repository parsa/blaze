//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/CompressedMatrix.h
//  \brief Specialization of the Creator class template for CompressedMatrix
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_COMPRESSEDMATRIX_H_
#define _BLAZETEST_MATHTEST_CREATOR_COMPRESSEDMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/creator/Default.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for compressed \f$ M \times N \f$ matrices.
//
// This specialization of the Creator class template is able to create random compressed matrices.
*/
template< typename T  // Element type of the compressed matrix
        , bool SO >   // Storage order of the compressed matrix
class Creator< blaze::CompressedMatrix<T,SO> >
{
 public:
   //**Type definitions****************************************************************************
   typedef blaze::CompressedMatrix<T,SO>  Type;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t m, size_t n, size_t nonzeros,
                            const Creator<T>& elementCreator = Creator<T>() );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   const blaze::CompressedMatrix<T,SO> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;         //!< The number of rows of the compressed matrix.
   size_t n_;         //!< The number of columns of the compressed matrix.
   size_t nonzeros_;  //!< The number of non-zero elements in the compressed matrix.
   Creator<T> ec_;    //!< Creator for the elements of the compressed matrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the creator specialization for CompressedMatrix.
//
// \param elementCreator The creator for the elements of the compressed matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename T  // Element type of the compressed matrix
        , bool SO >   // Storage order of the compressed matrix
inline Creator< blaze::CompressedMatrix<T,SO> >::Creator( const Creator<T>& elementCreator )
   : m_( 3UL )              // The number of rows of the compressed matrix
   , n_( 3UL )              // The number of columns of the compressed matrix
   , nonzeros_( 3UL )       // The total number of non-zero elements in the compressed matrix
   , ec_( elementCreator )  // Creator for the elements of the compressed matrix
{
   if( m_ * n_ < nonzeros_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for CompressedMatrix.
//
// \param m The number of rows of the compressed matrix.
// \param n The number of columns of the compressed matrix.
// \param nonzeros The number of non-zero elements in the compressed matrix.
// \param elementCreator The creator for the elements of the compressed matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename T  // Element type of the compressed matrix
        , bool SO >   // Storage order of the compressed matrix
inline Creator< blaze::CompressedMatrix<T,SO> >::Creator( size_t m, size_t n, size_t nonzeros,
                                                          const Creator<T>& elementCreator )
   : m_( m )                // The number of rows of the compressed matrix
   , n_( n )                // The number of columns of the compressed matrix
   , nonzeros_( nonzeros )  // The total number of non-zero elements in the compressed matrix
   , ec_( elementCreator )  // Creator for the elements of the compressed matrix
{
   if( m_ * n_ < nonzeros_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created compressed matrix.
//
// \return The randomly generated compressed matrix.
*/
template< typename T  // Element type of the compressed matrix
        , bool SO >   // Storage order of the compressed matrix
inline const blaze::CompressedMatrix<T,SO> Creator< blaze::CompressedMatrix<T,SO> >::operator()() const
{
   blaze::CompressedMatrix<T,SO> matrix( m_, n_, nonzeros_ );
   while( matrix.nonZeros() < nonzeros_ )
      matrix( blaze::rand<size_t>(0,m_-1), blaze::rand<size_t>(0,n_-1) ) = ec_();
   return matrix;
}
//*************************************************************************************************

} // namespace blazetest

#endif
