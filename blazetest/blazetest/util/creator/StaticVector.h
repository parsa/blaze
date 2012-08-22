//=================================================================================================
/*!
//  \file blazetest/util/creator/StaticVector.h
//  \brief Specialization of the Creator class template for StaticVector
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

#ifndef _BLAZETEST_UTIL_CREATOR_STATICVECTOR_H_
#define _BLAZETEST_UTIL_CREATOR_STATICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/StaticVector.h>
#include <blazetest/util/creator/Default.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for static vectors.
//
// This specialization of the Creator class template is able to create random static vectors.
*/
template< typename T  // Element type of the static vector
        , size_t N    // Number of elements of the static vector
        , bool TF >   // Transpose flag of the static vector
class Creator< blaze::StaticVector<T,N,TF> >
{
 public:
   //**Type definitions****************************************************************************
   typedef blaze::StaticVector<T,N,TF>  Type;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
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
   const blaze::StaticVector<T,N,TF> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Creator<T> ec_;  //!< Creator for the elements of the static vector.
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
/*!\brief Constructor for the creator specialization for StaticVector.
//
// \param elementCreator The creator for the elements of the static vector.
*/
template< typename T  // Element type of the static vector
        , size_t N    // Number of elements of the static vector
        , bool TF >   // Transpose flag of the static vector
inline Creator< blaze::StaticVector<T,N,TF> >::Creator( const Creator<T>& elementCreator )
   : ec_( elementCreator )  // Creator for the elements of the static vector
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created static vector.
//
// \return The randomly generated static vector.
*/
template< typename T  // Element type of the static vector
        , size_t N    // Number of elements of the static vector
        , bool TF >   // Transpose flag of the static vector
inline const blaze::StaticVector<T,N,TF> Creator< blaze::StaticVector<T,N,TF> >::operator()() const
{
   blaze::StaticVector<T,N,TF> vector;
   for( size_t i=0; i<N; ++i )
      vector[i] = ec_();
   return vector;
}
//*************************************************************************************************

} // namespace blazetest

#endif
