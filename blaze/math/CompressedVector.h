//=================================================================================================
/*!
//  \file blaze/math/CompressedVector.h
//  \brief Implementation of an arbitrarily sized compressed vector
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

#ifndef _BLAZE_MATH_COMPRESSEDVECTOR_H_
#define _BLAZE_MATH_COMPRESSEDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/sparse/CompressedMatrix.h>
#include <blaze/math/sparse/CompressedVector.h>
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
/*!\brief Specialization of the Rand class template for CompressedVector.
// \ingroup random
//
// This specialization of the Rand class creates random instances of CompressedVector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
class Rand< CompressedVector<Type,TF> >
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand( size_t size );
   explicit inline Rand( size_t size, size_t nonzeros );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator CompressedVector<Type,TF>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CompressedVector<Type,TF> vector_;  //!< The random vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Single argument constructor of the Rand specialization for CompressedVector.
//
// \param size The size of the random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline Rand< CompressedVector<Type,TF> >::Rand( size_t size )
   : vector_( size )  // The random vector
{
   if( size == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

   vector_.reserve( nonzeros );

   while( vector_.nonZeros() < nonzeros ) {
      vector_[ rand<size_t>( 0UL, size-1UL ) ] = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Two argument constructor of the Rand specialization for CompressedVector.
//
// \param size The size of the random vector.
// \param nonzeros The number of non-zero elements of the random vector.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline Rand< CompressedVector<Type,TF> >::Rand( size_t size, size_t nonzeros )
   : vector_( size, nonzeros )  // The random vector
{
   if( nonzeros > size )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( size == 0UL ) return;

   while( vector_.nonZeros() < nonzeros ) {
      vector_[ rand<size_t>( 0UL, size-1UL ) ] = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion to the created random CompressedVector.
//
// \return The random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline Rand< CompressedVector<Type,TF> >::operator CompressedVector<Type,TF>() const
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
/*!\brief Compressed single precision vector.
// \ingroup compressed_vector
*/
typedef CompressedVector<float,false>  CVecNf;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compressed double precision vector.
// \ingroup compressed_vector
*/
typedef CompressedVector<double,false>  CVecNd;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compressed vector with system-specific precision.
// \ingroup compressed_vector
*/
typedef CompressedVector<real,false>  CVecN;
//*************************************************************************************************

} // namespace blaze

#endif
