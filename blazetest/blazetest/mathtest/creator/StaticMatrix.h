//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/StaticMatrix.h
//  \brief Specialization of the Creator class template for StaticMatrix
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_STATICMATRIX_H_
#define _BLAZETEST_MATHTEST_CREATOR_STATICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/StaticMatrix.h>
#include <blazetest/mathtest/creator/Default.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for static matrices.
//
// This specialization of the Creator class template is able to create random static matrices.
*/
template< typename T  // Element type of the static matrix
        , size_t M    // Number of rows of the static matrix
        , size_t N    // Number of columns of the static matrix
        , bool SO >   // Storage order of the static matrix
class Creator< blaze::StaticMatrix<T,M,N,SO> >
{
 public:
   //**Type definitions****************************************************************************
   typedef blaze::StaticMatrix<T,M,N,SO>  Type;  //!< Type to be created by the Creator.
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
   const blaze::StaticMatrix<T,M,N,SO> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Creator<T> ec_;  //!< Creator for the elements of the static matrix.
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
/*!\brief Constructor for the creator specialization for StaticMatrix.
//
// \param elementCreator The creator for the elements of the static matrix.
*/
template< typename T  // Element type of the static matrix
        , size_t M    // Number of rows of the static matrix
        , size_t N    // Number of columns of the static matrix
        , bool SO >   // Storage order of the static matrix
inline Creator< blaze::StaticMatrix<T,M,N,SO> >::Creator( const Creator<T>& elementCreator )
   : ec_( elementCreator )  // Creator for the elements of the static matrix
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created static matrix.
//
// \return The randomly generated static matrix.
*/
template< typename T  // Element type of the static matrix
        , size_t M    // Number of rows of the static matrix
        , size_t N    // Number of columns of the static matrix
        , bool SO >   // Storage order of the static matrix
inline const blaze::StaticMatrix<T,M,N,SO> Creator< blaze::StaticMatrix<T,M,N,SO> >::operator()() const
{
   blaze::StaticMatrix<T,M,N,SO> matrix;

   // Initialization of a column-major matrix
   if( SO ) {
      for( size_t j=0UL; j<N; ++j )
         for( size_t i=0UL; i<M; ++i )
            matrix(i,j) = ec_();
   }

   // Initialization of a row-major matrix
   else {
      for( size_t i=0UL; i<M; ++i )
         for( size_t j=0UL; j<N; ++j )
            matrix(i,j) = ec_();
   }

   return matrix;
}
//*************************************************************************************************

} // namespace blazetest

#endif
