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
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline const DynamicVector<Type,TF> generate( size_t n ) const;
   inline const DynamicVector<Type,TF> generate( size_t n, Type min, Type max ) const;

   inline void randomize( DynamicVector<Type,TF>& vector ) const;
   inline void randomize( DynamicVector<Type,TF>& vector, Type min, Type max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random DynamicVector.
//
// \param n The size of the random vector.
// \return The generated random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const DynamicVector<Type,TF> Rand< DynamicVector<Type,TF> >::generate( size_t n ) const
{
   DynamicVector<Type,TF> vector( n );
   randomize( vector );
   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random DynamicVector.
//
// \param n The size of the random vector.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return The generated random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const DynamicVector<Type,TF> Rand< DynamicVector<Type,TF> >::generate( size_t n, Type min, Type max ) const
{
   DynamicVector<Type,TF> vector( n );
   randomize( vector, min, max );
   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a DynamicVector.
//
// \param vector The vector to be randomized.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void Rand< DynamicVector<Type,TF> >::randomize( DynamicVector<Type,TF>& vector ) const
{
   using blaze::randomize;

   const size_t size( vector.size() );
   for( size_t i=0UL; i<size; ++i ) {
      randomize( vector[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a DynamicVector.
//
// \param vector The vector to be randomized.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void Rand< DynamicVector<Type,TF> >::randomize( DynamicVector<Type,TF>& vector, Type min, Type max ) const
{
   using blaze::randomize;

   const size_t size( vector.size() );
   for( size_t i=0UL; i<size; ++i ) {
      randomize( vector[i], min, max );
   }
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
