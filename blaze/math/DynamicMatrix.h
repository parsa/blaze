//=================================================================================================
/*!
//  \file blaze/math/DynamicMatrix.h
//  \brief Header file for the complete DynamicMatrix implementation
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

#ifndef _BLAZE_MATH_DYNAMICMATRIX_H_
#define _BLAZE_MATH_DYNAMICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>
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
/*!\brief Specialization of the Rand class template for DynamicMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of DynamicMatrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
class Rand< DynamicMatrix<Type,SO> >
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand( size_t m, size_t n );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator DynamicMatrix<Type,SO>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   DynamicMatrix<Type,SO> matrix_;  //!< The random matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor of the Rand specialization for DynamicMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< DynamicMatrix<Type,SO> >::Rand( size_t m, size_t n )
   : matrix_( m, n )  // The random matrix
{
   for( size_t i=0UL; i<m; ++i ) {
      for( size_t j=0UL; j<n; ++j ) {
         matrix_(i,j) = rand<Type>();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion to the created random DynamicMatrix.
//
// \return The random matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Rand< DynamicMatrix<Type,SO> >::operator DynamicMatrix<Type,SO>() const
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
// \ingroup dynamic_matrix
*/
typedef DynamicMatrix<float,false>  MatMxNf;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief MxN double precision matrix.
// \ingroup dynamic_matrix
*/
typedef DynamicMatrix<double,false>  MatMxNd;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief MxN matrix with system-specific precision.
// \ingroup dynamic_matrix
*/
typedef DynamicMatrix<real,false>  MatMxN;
//*************************************************************************************************

} // namespace blaze

#endif
