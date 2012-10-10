//=================================================================================================
/*!
//  \file blaze/math/StaticMatrix.h
//  \brief Header file for the complete StaticMatrix implementation
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

#ifndef _BLAZE_MATH_STATICMATRIX_H_
#define _BLAZE_MATH_STATICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>
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
/*!\brief Specialization of the Rand class template for StaticMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of StaticMatrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
class Rand< StaticMatrix<Type,M,N,SO> >
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator StaticMatrix<Type,M,N,SO>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   StaticMatrix<Type,M,N,SO> matrix_;  //!< The random matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor of the Rand specialization for StaticMatrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline Rand< StaticMatrix<Type,M,N,SO> >::Rand()
   : matrix_()  // The random matrix
{
   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         matrix_(i,j) = rand<Type>();
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion to the created random StaticMatrix.
//
// \return The random matrix.
*/
template< typename Type  // Data type of the matrix
        , size_t M       // Number of rows
        , size_t N       // Number of columns
        , bool SO >      // Storage order
inline Rand< StaticMatrix<Type,M,N,SO> >::operator StaticMatrix<Type,M,N,SO>() const
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
/*!\brief 2x2 single precision matrix.
// \ingroup static_matrix_2x2
*/
typedef StaticMatrix<float,2UL,2UL,false>  Mat2x2f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2x2 double precision matrix.
// \ingroup static_matrix_2x2
*/
typedef StaticMatrix<double,2UL,2UL,false>  Mat2x2d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2x2 matrix with system-specific precision.
// \ingroup static_matrix_2x2
*/
typedef StaticMatrix<real,2UL,2UL,false>  Mat2x2;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 3x3 single precision matrix.
// \ingroup static_matrix_3x3
*/
typedef StaticMatrix<float,3UL,3UL,false>  Mat3x3f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 3x3 double precision matrix.
// \ingroup static_matrix_3x3
*/
typedef StaticMatrix<double,3UL,3UL,false>  Mat3x3d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 3x3 matrix with system-specific precision.
// \ingroup static_matrix_3x3
*/
typedef StaticMatrix<real,3UL,3UL,false>  Mat3x3;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 4x4 single precision matrix.
// \ingroup static_matrix_4x4
*/
typedef StaticMatrix<float,4UL,4UL,false>  Mat4x4f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 4x4 double precision matrix.
// \ingroup static_matrix_4x4
*/
typedef StaticMatrix<double,4UL,4UL,false>  Mat4x4d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 4x4 matrix with system-specific precision.
// \ingroup static_matrix_4x4
*/
typedef StaticMatrix<real,4UL,4UL,false>  Mat4x4;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 5x5 single precision matrix.
// \ingroup static_matrix_5x5
*/
typedef StaticMatrix<float,5UL,5UL,false>  Mat5x5f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 5x5 double precision matrix.
// \ingroup static_matrix_5x5
*/
typedef StaticMatrix<double,5UL,5UL,false>  Mat5x5d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 5x5 matrix with system-specific precision.
// \ingroup static_matrix_5x5
*/
typedef StaticMatrix<real,5UL,5UL,false>  Mat5x5;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 6x6 single precision matrix.
// \ingroup static_matrix_6x6
*/
typedef StaticMatrix<float,6UL,6UL,false>  Mat6x6f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 6x6 double precision matrix.
// \ingroup static_matrix_6x6
*/
typedef StaticMatrix<double,6UL,6UL,false>  Mat6x6d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 6x6 matrix with system-specific precision.
// \ingroup static_matrix_6x6
*/
typedef StaticMatrix<real,6UL,6UL,false>  Mat6x6;
//*************************************************************************************************

} // namespace blaze

#endif
