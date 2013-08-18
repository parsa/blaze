//=================================================================================================
/*!
//  \file blaze/math/StaticVector.h
//  \brief Header file for the complete StaticVector implementation
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

#ifndef _BLAZE_MATH_STATICVECTOR_H_
#define _BLAZE_MATH_STATICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/DenseVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticMatrix.h>
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
/*!\brief Specialization of the Rand class template for StaticVector.
// \ingroup random
//
// This specialization of the Rand class creates random instances of StaticVector.
*/
template< typename Type  // Data type of the vector
        , size_t N       // Number of elements
        , bool TF >      // Transpose flag
class Rand< StaticVector<Type,N,TF> >
{
 public:
   //**Generate functions**************************************************************************
   /*!\name Generate functions */
   //@{
   inline const StaticVector<Type,N,TF> generate() const;

   template< typename Arg >
   inline const StaticVector<Type,N,TF> generate( const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************

   //**Randomize functions*************************************************************************
   /*!\name Randomize functions */
   //@{
   inline void randomize( StaticVector<Type,N,TF>& vector ) const;

   template< typename Arg >
   inline void randomize( StaticVector<Type,N,TF>& vector, const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random StaticVector.
//
// \return The generated random vector.
*/
template< typename Type  // Data type of the vector
        , size_t N       // Number of elements
        , bool TF >      // Transpose flag
inline const StaticVector<Type,N,TF> Rand< StaticVector<Type,N,TF> >::generate() const
{
   StaticVector<Type,N,TF> vector;
   randomize( vector );
   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random StaticVector.
//
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return The generated random vector.
*/
template< typename Type   // Data type of the vector
        , size_t N        // Number of elements
        , bool TF >       // Transpose flag
template< typename Arg >  // Min/max argument type
inline const StaticVector<Type,N,TF>
   Rand< StaticVector<Type,N,TF> >::generate( const Arg& min, const Arg& max ) const
{
   StaticVector<Type,N,TF> vector;
   randomize( vector, min, max );
   return vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a StaticVector.
//
// \param vector The vector to be randomized.
// \return void
*/
template< typename Type  // Data type of the vector
        , size_t N       // Number of elements
        , bool TF >      // Transpose flag
inline void Rand< StaticVector<Type,N,TF> >::randomize( StaticVector<Type,N,TF>& vector ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<N; ++i ) {
      randomize( vector[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a StaticVector.
//
// \param vector The vector to be randomized.
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
// \return void
*/
template< typename Type   // Data type of the vector
        , size_t N        // Number of elements
        , bool TF >       // Transpose flag
template< typename Arg >  // Min/max argument type
inline void Rand< StaticVector<Type,N,TF> >::randomize( StaticVector<Type,N,TF>& vector,
                                                        const Arg& min, const Arg& max ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<N; ++i ) {
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
/*!\brief 2-dimensional single precision vector.
// \ingroup static_vector_2
*/
typedef StaticVector<float,2UL,false>  Vec2f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2-dimensional double precision vector.
// \ingroup static_vector_2
*/
typedef StaticVector<double,2UL,false>  Vec2d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2-dimensional vector with system-specific precision.
// \ingroup static_vector_2
*/
typedef StaticVector<real,2UL,false>  Vec2;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 3-dimensional single precision vector.
// \ingroup static_vector_3
*/
typedef StaticVector<float,3UL,false>  Vec3f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 3-dimensional double precision vector.
// \ingroup static_vector_3
*/
typedef StaticVector<double,3UL,false>  Vec3d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 3-dimensional vector with system-specific precision.
// \ingroup static_vector_3
*/
typedef StaticVector<real,3UL,false>  Vec3;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 6-dimensional single precision vector.
// \ingroup static_vector_6
*/
typedef StaticVector<float,6UL,false>  Vec6f;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 6-dimensional double precision vector.
// \ingroup static_vector_6
*/
typedef StaticVector<double,6UL,false>  Vec6d;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 6-dimensional vector with system-specific precision.
// \ingroup static_vector_6
*/
typedef StaticVector<real,6UL,false>  Vec6;
//*************************************************************************************************

} // namespace blaze

#endif
