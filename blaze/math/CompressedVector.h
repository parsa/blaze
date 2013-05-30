//=================================================================================================
/*!
//  \file blaze/math/CompressedVector.h
//  \brief Header file for the complete CompressedVector implementation
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

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/sparse/CompressedVector.h>
#include <blaze/math/SparseVector.h>
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
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline const CompressedVector<Type,TF> generate( size_t size ) const;
   inline const CompressedVector<Type,TF> generate( size_t size, size_t nonzeros ) const;
   inline const CompressedVector<Type,TF> generate( size_t size, Type min, Type max ) const;
   inline const CompressedVector<Type,TF> generate( size_t size, size_t nonzeros, Type min, Type max ) const;

   inline void randomize( CompressedVector<Type,TF>& vector ) const;
   inline void randomize( CompressedVector<Type,TF>& vector, size_t nonzeros ) const;
   inline void randomize( CompressedVector<Type,TF>& vector, Type min, Type max ) const;
   inline void randomize( CompressedVector<Type,TF>& vector, size_t nonzeros, Type min, Type max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedVector.
//
// \param size The size of the random vector.
// \return The generated random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const CompressedVector<Type,TF>
   Rand< CompressedVector<Type,TF> >::generate( size_t size ) const
{
   if( size == 0UL ) return;

   CompressedVector<Type,TF> vector( size );
   randomize( vector );

   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedVector.
//
// \param size The size of the random vector.
// \param nonzeros The number of non-zero elements of the random vector.
// \return The generated random vector.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const CompressedVector<Type,TF>
   Rand< CompressedVector<Type,TF> >::generate( size_t size, size_t nonzeros ) const
{
   if( nonzeros > size )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( size == 0UL ) return;

   CompressedVector<Type,TF> vector( size, nonzeros );
   randomize( vector, nonzeros );

   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedVector.
//
// \param size The size of the random vector.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return The generated random vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const CompressedVector<Type,TF>
   Rand< CompressedVector<Type,TF> >::generate( size_t size, Type min, Type max ) const
{
   if( size == 0UL ) return;

   CompressedVector<Type,TF> vector( size );
   randomize( vector, min, max );

   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedVector.
//
// \param size The size of the random vector.
// \param nonzeros The number of non-zero elements of the random vector.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return The generated random vector.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const CompressedVector<Type,TF>
   Rand< CompressedVector<Type,TF> >::generate( size_t size, size_t nonzeros, Type min, Type max ) const
{
   if( nonzeros > size )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( size == 0UL ) return;

   CompressedVector<Type,TF> vector( size, nonzeros );
   randomize( vector, nonzeros, min, max );

   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedVector.
//
// \param vector The vector to be randomized.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void Rand< CompressedVector<Type,TF> >::randomize( CompressedVector<Type,TF>& vector ) const
{
   const size_t size( vector.size() );

   if( size == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

   vector.reset();
   vector.reserve( nonzeros );

   while( vector.nonZeros() < nonzeros ) {
      vector[ rand<size_t>( 0UL, size-1UL ) ] = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedVector.
//
// \param vector The vector to be randomized.
// \param nonzeros The number of non-zero elements of the random vector.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void Rand< CompressedVector<Type,TF> >::randomize( CompressedVector<Type,TF>& vector, size_t nonzeros ) const
{
   const size_t size( vector.size() );

   if( nonzeros > size )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( size == 0UL ) return;

   vector.reset();
   vector.reserve( nonzeros );

   while( vector.nonZeros() < nonzeros ) {
      vector[ rand<size_t>( 0UL, size-1UL ) ] = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedVector.
//
// \param vector The vector to be randomized.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void Rand< CompressedVector<Type,TF> >::randomize( CompressedVector<Type,TF>& vector, Type min, Type max ) const
{
   const size_t size( vector.size() );

   if( size == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

   vector.reset();
   vector.reserve( nonzeros );

   while( vector.nonZeros() < nonzeros ) {
      vector[ rand<size_t>( 0UL, size-1UL ) ] = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedVector.
//
// \param vector The vector to be randomized.
// \param nonzeros The number of non-zero elements of the random vector.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void Rand< CompressedVector<Type,TF> >::randomize( CompressedVector<Type,TF>& vector, size_t nonzeros, Type min, Type max ) const
{
   const size_t size( vector.size() );

   if( nonzeros > size )
      throw std::invalid_argument( "Invalid number of non-zero elements" );

   if( size == 0UL ) return;

   vector.reset();
   vector.reserve( nonzeros );

   while( vector.nonZeros() < nonzeros ) {
      vector[ rand<size_t>( 0UL, size-1UL ) ] = rand<Type>( min, max );
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
