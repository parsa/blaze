//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/CompressedVector.h
//  \brief Specialization of the Creator class template for CompressedVector
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_COMPRESSEDVECTOR_H_
#define _BLAZETEST_MATHTEST_CREATOR_COMPRESSEDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/CompressedVector.h>
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
/*!\brief Specialization of the Creator class template for N-dimensional compressed vectors.
//
// This specialization of the Creator class template is able to create random N-dimensional
// compressed vectors.
*/
template< typename T  // Element type of the N-dimensional compressed vector
        , bool TF >   // Transpose flag of the N-dimensional compressed vector
class Creator< blaze::CompressedVector<T,TF> >
{
 public:
   //**Type definitions****************************************************************************
   typedef blaze::CompressedVector<T,TF>  Type;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t size, size_t nonzeros,
                            const Creator<T>& elementCreator = Creator<T>() );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   const blaze::CompressedVector<T,TF> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;      //!< The size for the N-dimensional compressed vector.
   size_t nonzeros_;  //!< The number of non-zero elements in the compressed vector.
   Creator<T> ec_;    //!< Creator for the elements of the N-dimensional compressed vector.
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
/*!\brief Constructor for the creator specialization for CompressedVector.
//
// \param elementCreator The creator for the elements of the N-dimensional compressed vector.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename T  // Element type of the N-dimensional compressed vector
        , bool TF >   // Transpose flag of the N-dimensional compressed vector
inline Creator< blaze::CompressedVector<T,TF> >::Creator( const Creator<T>& elementCreator )
   : size_( 3UL )           // The size for the N-dimensional compressed vector
   , nonzeros_( 1UL )       // The number of non-zero elements in the compressed vector
   , ec_( elementCreator )  // Creator for the elements of the N-dimensional compressed vector
{
   if( size_ < nonzeros_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for CompressedVector.
//
// \param size The size for the N-dimensional compressed vector.
// \param nonzeros The number of non-zero elements in the compressed vector.
// \param elementCreator The creator for the elements of the N-dimensional compressed vector.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename T  // Element type of the N-dimensional compressed vector
        , bool TF >   // Transpose flag of the N-dimensional compressed vector
inline Creator< blaze::CompressedVector<T,TF> >::Creator( size_t size, size_t nonzeros, const Creator<T>& elementCreator )
   : size_( size )          // The size for the N-dimensional compressed vector
   , nonzeros_( nonzeros )  // The number of non-zero elements in the compressed vector
   , ec_( elementCreator )  // Creator for the elements of the N-dimensional compressed vector
{
   if( size_ < nonzeros_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created N-dimensional compressed vector.
//
// \return The randomly generated N-dimensional compressed vector.
*/
template< typename T  // Element type of the N-dimensional compressed vector
        , bool TF >   // Transpose flag of the N-dimensional compressed vector
inline const blaze::CompressedVector<T,TF> Creator< blaze::CompressedVector<T,TF> >::operator()() const
{
   blaze::CompressedVector<T,TF> vector( size_ );
   while( vector.nonZeros() < nonzeros_ )
      vector[blaze::rand<size_t>(0,size_-1)] = ec_();
   return vector;
}
//*************************************************************************************************

} // namespace blazetest

#endif
