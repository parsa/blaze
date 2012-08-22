//=================================================================================================
/*!
//  \file blazetest/util/creator/DynamicVector.h
//  \brief Specialization of the Creator class template for DynamicVector
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

#ifndef _BLAZETEST_UTIL_CREATOR_DYNAMICVECTOR_H_
#define _BLAZETEST_UTIL_CREATOR_DYNAMICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/DynamicVector.h>
#include <blazetest/util/creator/Default.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for N-dimensional vectors.
//
// This specialization of the Creator class template is able to create random N-dimensional
// vectors.
*/
template< typename T  // Element type of the N-dimensional vector
        , bool TF >   // Transpose flag of the N-dimensional vector
class Creator< blaze::DynamicVector<T,TF> >
{
 public:
   //**Type definitions****************************************************************************
   typedef blaze::DynamicVector<T,TF>  Type;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t size, const Creator<T>& elementCreator = Creator<T>() );
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
   const blaze::DynamicVector<T,TF> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;    //!< The size for the N-dimensional vector.
   Creator<T> ec_;  //!< Creator for the elements of the N-dimensional vector.
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
/*!\brief Constructor for the creator specialization for DynamicVector.
//
// \param elementCreator The creator for the elements of the N-dimensional vector.
*/
template< typename T  // Element type of the N-dimensional vector
        , bool TF >   // Transpose flag of the N-dimensional vector
inline Creator< blaze::DynamicVector<T,TF> >::Creator( const Creator<T>& elementCreator )
   : size_( 3UL )           // The size for the N-dimensional vector
   , ec_( elementCreator )  // Creator for the elements of the N-dimensional vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for DynamicVector.
//
// \param size The size for the N-dimensional vector.
// \param elementCreator The creator for the elements of the N-dimensional vector.
*/
template< typename T  // Element type of the N-dimensional vector
        , bool TF >   // Transpose flag of the N-dimensional vector
inline Creator< blaze::DynamicVector<T,TF> >::Creator( size_t size, const Creator<T>& elementCreator )
   : size_( size )          // The size for the N-dimensional vector
   , ec_( elementCreator )  // Creator for the elements of the N-dimensional vector
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created N-dimensional vector.
//
// \return The randomly generated N-dimensional vector.
*/
template< typename T  // Element type of the N-dimensional vector
        , bool TF >   // Transpose flag of the N-dimensional vector
inline const blaze::DynamicVector<T,TF> Creator< blaze::DynamicVector<T,TF> >::operator()() const
{
   blaze::DynamicVector<T,TF> vector( size_ );
   for( size_t i=0; i<size_; ++i )
      vector[i] = ec_();
   return vector;
}
//*************************************************************************************************

} // namespace blazetest

#endif
