//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/DynamicMatrix.h
//  \brief Specialization of the Creator class template for DynamicMatrix
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_DYNAMICMATRIX_H_
#define _BLAZETEST_MATHTEST_CREATOR_DYNAMICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/DynamicMatrix.h>
#include <blazetest/mathtest/creator/Default.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for dynamic \f$ M \times N \f$ matrices.
//
// This specialization of the Creator class template is able to create random \f$ M \times N \f$
// matrices.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
class Creator< blaze::DynamicMatrix<T,SO> >
{
 public:
   //**Type definitions****************************************************************************
   typedef blaze::DynamicMatrix<T,SO>  Type;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t m, size_t n, const Creator<T>& elementCreator = Creator<T>() );
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
   const blaze::DynamicMatrix<T,SO> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;       //!< The number of rows of the dynamic matrix.
   size_t n_;       //!< The number of columns of the dynamic matrix.
   Creator<T> ec_;  //!< Creator for the elements of the dynamic matrix.
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
/*!\brief Constructor for the creator specialization for DynamicMatrix.
//
// \param elementCreator The creator for the elements of the dynamic matrix.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
inline Creator< blaze::DynamicMatrix<T,SO> >::Creator( const Creator<T>& elementCreator )
   : m_( 3UL )              // The number of rows of the dynamic matrix
   , n_( 3UL )              // The number of columns of the dynamic matrix
   , ec_( elementCreator )  // Creator for the elements of the dynamic matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for DynamicMatrix.
//
// \param m The number of rows of the dynamic matrix.
// \param n The number of columns of the dynamic matrix.
// \param elementCreator The creator for the elements of the dynamic matrix.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
inline Creator< blaze::DynamicMatrix<T,SO> >::Creator( size_t m, size_t n, const Creator<T>& elementCreator )
   : m_( m )                // The number of rows of the dynamic matrix
   , n_( n )                // The number of columns of the dynamic matrix
   , ec_( elementCreator )  // Creator for the elements of the dynamic matrix
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created dynamic matrix.
//
// \return The randomly generated dynamic matrix.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
inline const blaze::DynamicMatrix<T,SO> Creator< blaze::DynamicMatrix<T,SO> >::operator()() const
{
   blaze::DynamicMatrix<T,SO> matrix( m_, n_ );

   // Initialization of a column-major matrix
   if( SO ) {
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=0UL; i<m_; ++i )
            matrix(i,j) = ec_();
   }

   // Initialization of a row-major matrix
   else {
      for( size_t i=0UL; i<m_; ++i )
         for( size_t j=0UL; j<n_; ++j )
            matrix(i,j) = ec_();
   }

   return matrix;
}
//*************************************************************************************************

} // namespace blazetest

#endif
