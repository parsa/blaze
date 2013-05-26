//=================================================================================================
/*!
//  \file blaze/math/DynamicVector.h
//  \brief Header file for the complete DynamicVector implementation
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

#ifndef _BLAZE_MATH_DYNAMICVECTOR_H_
#define _BLAZE_MATH_DYNAMICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/DenseVector.h>
#include <blaze/math/DynamicMatrix.h>
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
/*!\brief Specialization of the Rand class template for DynamicVector.
// \ingroup random
//
// This specialization of the Rand class creates random instances of DynamicVector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
class Rand< DynamicVector<Type,TF> >
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand( size_t n );
   explicit inline Rand( size_t n, Type min, Type max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator DynamicVector<Type,TF>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   DynamicVector<Type,TF> vector_;  //!< The random vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor of the Rand specialization for DynamicVector.
//
// \param n The size of the random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline Rand< DynamicVector<Type,TF> >::Rand( size_t n )
   : vector_( n )  // The random vector
{
   for( size_t i=0UL; i<n; ++i ) {
      vector_[i] = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand specialization for DynamicVector.
//
// \param n The size of the random vector.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline Rand< DynamicVector<Type,TF> >::Rand( size_t n, Type min, Type max )
   : vector_( n )  // The random vector
{
   for( size_t i=0UL; i<n; ++i ) {
      vector_[i] = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random DynamicVector.
//
// \return The random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline Rand< DynamicVector<Type,TF> >::operator DynamicVector<Type,TF>() const
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief N-dimensional single precision vector.
// \ingroup dynamic_vector
*/
typedef DynamicVector<float,false>  VecNf;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief N-dimensional double precision vector.
// \ingroup dynamic_vector
*/
typedef DynamicVector<double,false>  VecNd;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief N-dimensional vector with system-specific precision.
// \ingroup dynamic_vector
*/
typedef DynamicVector<real,false>  VecN;
//*************************************************************************************************

} // namespace blaze

#endif
